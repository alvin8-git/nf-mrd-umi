#!/usr/bin/env python3
"""validate.py — offline validation harness for the MRD engine (CLSI EP17-flavored).

Breaks the circular-self-test trap. The Part-3 `selftest`s draw data under the
engine's OWN null and then confirm it passes — they prove code fidelity, not
model validity. This harness does the opposite on purpose:

  1. Cohort hygiene — fits the per-site background on a TRAINING fold of blanks
     and evaluates a DISJOINT held-out fold. You never validate specificity on
     the donors you trained the null on.
  2. Adversarial truth — simulates data under a RICHER model than the engine
     assumes: a per-sample shared batch multiplier that makes all sites move
     together (cross-site error correlation). The engine's independent
     Monte-Carlo null cannot see this, so the harness can expose the
     anticonservative-null failure as an inflated false-positive rate.

Outputs:
  - Limit of Blank (LoB) + observed false-positive rate vs nominal alpha  (specificity)
  - per-VAF hit-rate, LoD95, LoQ, quantification bias                     (sensitivity)
  - enrichment recovery (calibrated per-probe enrich vs truth)

simulate-run : end-to-end on simulated cohorts, no real data needed
run-real     : evaluate real per-sample site-count TSVs listed in a manifest
selftest     : assert harness mechanics AND that it detects null inflation
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mrd_integrate import (Site, call_mrd, fit_betabinom_mom,          # noqa: E402
                           fit_enrichment_from_dilution, load_sites,
                           load_donor_rates)


# ---------------------------------------------------------------------------
# Simulation model (deliberately richer than the engine's null)
# ---------------------------------------------------------------------------
@dataclass
class ProbeSpec:
    site_id: str
    error_rate: float   # true per-site background error rate
    enrich: float       # true per-probe allele-specific enrichment factor


def make_panel(rng, n: int, enrich: float, mean_err: float = 5e-4) -> list[ProbeSpec]:
    """Per-site error rates vary log-normally (real panels are heterogeneous)."""
    rates = np.exp(rng.normal(np.log(mean_err), 0.5, size=n))
    return [ProbeSpec(f"site{i}", float(rates[i]), enrich) for i in range(n)]


def simulate_counts(rng, panel, vaf: float, batch_sd: float,
                    depth: int) -> dict[str, tuple[int, int]]:
    """One sample's per-site (depth, alt). `batch_sd`>0 applies a per-sample
    shared multiplier to every site's error rate -> cross-site CORRELATION that
    the engine's independent null does not model. Tumor signal is enriched."""
    m = float(np.exp(rng.normal(0.0, batch_sd))) if batch_sd > 0 else 1.0
    out = {}
    for p in panel:
        err = int(rng.binomial(depth, min(0.5, p.error_rate * m)))
        if vaf > 0:                              # odds-model enrichment
            o = vaf / (1.0 - vaf)
            obs_p = (p.enrich * o) / (1.0 + p.enrich * o)
            tum = int(rng.binomial(depth, obs_p))
        else:
            tum = 0
        out[p.site_id] = (depth, min(depth, err + tum))
    return out


# ---------------------------------------------------------------------------
# Build the engine's inputs (background + enrichment) from TRAINING data only
# ---------------------------------------------------------------------------
def build_background(rng, panel, n_train: int, depth: int,
                     batch_sd: float) -> dict[str, tuple[float, float]]:
    """MoM beta-binomial per site over n_train BLANK training samples."""
    alt = {p.site_id: [] for p in panel}
    dep = {p.site_id: [] for p in panel}
    for _ in range(n_train):
        c = simulate_counts(rng, panel, 0.0, batch_sd, depth)
        for sid, (d, a) in c.items():
            alt[sid].append(a)
            dep[sid].append(d)
    bg, donor_rates = {}, {}
    for p in panel:
        av, dv = np.array(alt[p.site_id]), np.array(dep[p.site_id])
        bg[p.site_id] = fit_betabinom_mom(av, dv)
        # per-donor error-rate vector keyed to the Site key used downstream;
        # this matrix IS the empirical null (preserves cross-site covariance).
        donor_rates[(p.site_id, 0, "C", "T")] = av / np.maximum(dv, 1)
    return bg, donor_rates


def build_enrichment(rng, panel, bg, depth: int,
                     levels=(1e-3, 5e-3, 1e-2, 5e-2),
                     reps: int = 5) -> dict[str, float]:
    """Calibrate per-probe enrichment from a contrived dilution series, the same
    way the real assay would — exercises fit_enrichment_from_dilution. The engine
    then gets the CALIBRATED (not the true) enrichment, as in production."""
    enrich = {}
    for p in panel:
        bg_rate = bg[p.site_id][0] / (bg[p.site_id][0] + bg[p.site_id][1])
        known, observed = [], []
        for v in levels:
            for _ in range(reps):
                _, a = simulate_counts(rng, [p], v, 0.0, depth)[p.site_id]
                known.append(v)
                observed.append(max(0.0, (a - depth * bg_rate) / depth))
        enrich[p.site_id] = fit_enrichment_from_dilution(np.array(known),
                                                         np.array(observed))
    return enrich


def counts_to_sites(counts, bg, enrich) -> list[Site]:
    sites = []
    for sid, (d, a) in counts.items():
        alpha, beta = bg[sid]
        sites.append(Site(sid, 0, "C", "T", d, a, alpha, beta, enrich[sid]))
    return sites


# ---------------------------------------------------------------------------
# Metrics (CLSI EP17-flavored)
# ---------------------------------------------------------------------------
def eval_samples(sites_list, rng, alpha, min_molecules, n_iter, n_boot,
                 donor_rates=None):
    calls, tfs = [], []
    for sites in sites_list:
        r = call_mrd(sites, rng, p_threshold=alpha,
                     min_total_molecules=min_molecules,
                     n_iter=n_iter, n_boot=n_boot, donor_rates=donor_rates)
        calls.append(r["call"])
        tfs.append(r["tumor_fraction"])
    return calls, tfs


def lob_and_fpr(blank_calls, blank_tfs):
    pos = sum(c == "positive" for c in blank_calls)
    fpr = pos / len(blank_calls) if blank_calls else float("nan")
    lob = float(np.percentile(blank_tfs, 95)) if blank_tfs else float("nan")
    return fpr, lob


def lod_loq(level_results, bias_tol=0.3, cv_tol=0.3):
    """level_results: list of (vaf, calls, tfs). Returns per-level table, LoD95, LoQ."""
    table, lod95, loq = [], None, None
    for vaf, calls, tfs in sorted(level_results, key=lambda x: x[0]):
        hit = sum(c == "positive" for c in calls) / len(calls)
        mean_tf = float(np.mean(tfs))
        cv = float(np.std(tfs) / mean_tf) if mean_tf > 0 else float("inf")
        rel_bias = abs(mean_tf - vaf) / vaf if vaf > 0 else float("inf")
        table.append({"vaf": vaf, "hit_rate": round(hit, 3),
                      "mean_tf": mean_tf, "cv": round(cv, 3),
                      "rel_bias": round(rel_bias, 3)})
        if lod95 is None and hit >= 0.95:
            lod95 = vaf
        if loq is None and rel_bias <= bias_tol and cv <= cv_tol:
            loq = vaf
    return table, lod95, loq


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------
def run_simulation(cfg, rng):
    panel = make_panel(rng, cfg["panel_size"], cfg["enrich"])
    bg, donor_rates = build_background(rng, panel, cfg["n_train"], cfg["depth"],
                                       cfg["batch_sd"])
    enrich = build_enrichment(rng, panel, bg, cfg["depth"])
    # Empirical (covariance-preserving) null by default; None = broken independent.
    dr = donor_rates if cfg.get("use_empirical", True) else None

    # Held-out blanks (disjoint from training) -> specificity
    blank_sites = [counts_to_sites(simulate_counts(rng, panel, 0.0,
                   cfg["batch_sd"], cfg["depth"]), bg, enrich)
                   for _ in range(cfg["n_blank"])]
    bc, btf = eval_samples(blank_sites, rng, cfg["alpha"], cfg["min_molecules"],
                           cfg["n_iter"], cfg["n_boot"], donor_rates=dr)
    fpr, lob = lob_and_fpr(bc, btf)

    # Dilution series -> sensitivity
    level_results = []
    for v in cfg["vaf_levels"]:
        sites_l = [counts_to_sites(simulate_counts(rng, panel, v, cfg["batch_sd"],
                   cfg["depth"]), bg, enrich) for _ in range(cfg["n_rep"])]
        calls, tfs = eval_samples(sites_l, rng, cfg["alpha"], cfg["min_molecules"],
                                  cfg["n_iter"], cfg["n_boot"], donor_rates=dr)
        level_results.append((v, calls, tfs))
    table, lod95, loq = lod_loq(level_results)

    enrich_recovery = float(np.mean([enrich[p.site_id] for p in panel])) / cfg["enrich"]
    return {
        "config": cfg,
        "null_model": "empirical_donor_resampling" if dr is not None
                      else "independent_betabinom",
        "nominal_alpha": cfg["alpha"],
        "observed_fpr": round(fpr, 4),
        "specificity": round(1 - fpr, 4),
        "limit_of_blank_tf": lob,
        "lod95_vaf": lod95,
        "loq_vaf": loq,
        "dilution": table,
        "enrichment_recovery_ratio": round(enrich_recovery, 3),
        "specificity_holds": bool(fpr <= cfg["alpha"] + 0.02),
    }


def write_report(report, path):
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2)


def print_summary(r):
    print(f"  null model        : {r.get('null_model', '?')}")
    print(f"  nominal alpha     : {r['nominal_alpha']}")
    print(f"  observed FPR      : {r['observed_fpr']}   "
          f"({'OK' if r['specificity_holds'] else 'INFLATED — null anticonservative'})")
    print(f"  LoB (tumor frac)  : {r['limit_of_blank_tf']:.2e}")
    print(f"  LoD95 VAF         : {r['lod95_vaf']}")
    print(f"  LoQ VAF           : {r['loq_vaf']}")
    print(f"  enrichment recov. : {r['enrichment_recovery_ratio']}x of truth")
    for row in r["dilution"]:
        print(f"    VAF {row['vaf']:.4f}  hit={row['hit_rate']:.2f}  "
              f"meanTF={row['mean_tf']:.2e}  CV={row['cv']:.2f}  bias={row['rel_bias']:.2f}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _cfg_from_args(args):
    return {
        "panel_size": args.panel_size, "enrich": args.enrich,
        "depth": args.depth, "batch_sd": args.batch_sd,
        "n_train": args.n_train, "n_blank": args.n_blank, "n_rep": args.n_rep,
        "vaf_levels": args.vaf_levels, "alpha": args.alpha,
        "min_molecules": args.min_molecules,
        "n_iter": args.iterations, "n_boot": args.bootstrap,
        "use_empirical": args.null == "empirical",
    }


def cmd_simulate(args):
    rng = np.random.default_rng(args.seed)
    r = run_simulation(_cfg_from_args(args), rng)
    write_report(r, args.out)
    print(f"validation report -> {args.out}")
    print_summary(r)


def cmd_run_real(args):
    """Evaluate real samples. Manifest TSV columns: sample, label[blank|tumor],
    vaf_truth(optional), site_counts_path. Shared --background and the engine's
    enrichment column score every sample the production way."""
    rng = np.random.default_rng(args.seed)
    donor_rates = load_donor_rates(args.pon) if args.pon else None
    blanks, levels = [], {}
    with open(args.manifest, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sites = load_sites(row["site_counts_path"], args.background)
            r = call_mrd(sites, rng, p_threshold=args.alpha,
                         min_total_molecules=args.min_molecules,
                         n_iter=args.iterations, n_boot=args.bootstrap,
                         donor_rates=donor_rates)
            if row["label"] == "blank":
                blanks.append(r)
            else:
                v = float(row.get("vaf_truth") or 0.0)
                levels.setdefault(v, []).append(r)
    fpr, lob = lob_and_fpr([b["call"] for b in blanks],
                           [b["tumor_fraction"] for b in blanks])
    level_results = [(v, [r["call"] for r in rs], [r["tumor_fraction"] for r in rs])
                     for v, rs in levels.items()]
    table, lod95, loq = lod_loq(level_results) if level_results else ([], None, None)
    r = {"nominal_alpha": args.alpha, "observed_fpr": round(fpr, 4),
         "specificity": round(1 - fpr, 4), "limit_of_blank_tf": lob,
         "lod95_vaf": lod95, "loq_vaf": loq, "dilution": table,
         "n_blanks": len(blanks), "specificity_holds": bool(fpr <= args.alpha + 0.02)}
    write_report(r, args.out)
    print(f"validation report -> {args.out}")
    print_summary({**r, "enrichment_recovery_ratio": float("nan")})


def cmd_selftest(_args):
    base = dict(panel_size=20, enrich=30.0, depth=4000, n_train=120, n_blank=120,
                n_rep=12, vaf_levels=[1e-3, 3e-3, 1e-2], alpha=0.05,
                min_molecules=1000, n_iter=1500, n_boot=40, batch_sd=0.7)

    # SAME correlated-batch data, two nulls: independent (broken) vs empirical (fix).
    r_indep = run_simulation({**base, "use_empirical": False}, np.random.default_rng(7))
    r_emp = run_simulation({**base, "use_empirical": True}, np.random.default_rng(7))

    # Same data, two nulls. Empirical (fatter tail) is structurally no worse, and
    # brings the correlated-batch FPR back near nominal; independent exceeds it.
    assert r_indep["observed_fpr"] > base["alpha"], r_indep   # independent inflates
    assert r_emp["observed_fpr"] <= r_indep["observed_fpr"], \
        (r_indep["observed_fpr"], r_emp["observed_fpr"])
    assert r_emp["observed_fpr"] <= 0.12, r_emp              # empirical controlled

    # The fix must not cost sensitivity: strong signal still detected.
    top = [row for row in r_emp["dilution"] if row["vaf"] == 1e-2][0]
    assert top["hit_rate"] >= 0.9, top
    # Odds-space enrichment calibration is now ~unbiased -> quantification valid.
    assert 0.8 <= r_emp["enrichment_recovery_ratio"] <= 1.2, \
        r_emp["enrichment_recovery_ratio"]
    assert top["rel_bias"] <= 0.3, top                # TF tracks truth now
    assert r_emp["loq_vaf"] is not None, r_emp        # LoQ exists (was None)

    print("validate.py selftest: PASS")
    print(f"  independent null FPR = {r_indep['observed_fpr']} "
          f"(broken; nominal {base['alpha']})")
    print(f"  empirical   null FPR = {r_emp['observed_fpr']} "
          f"(fixed; covariance-preserving)")
    print(f"  enrichment recovery = {r_emp['enrichment_recovery_ratio']}x  "
          f"LoQ = {r_emp['loq_vaf']}  (was None)")
    print(f"  empirical LoD95 = {r_emp['lod95_vaf']}  hit@1% = {top['hit_rate']}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate-run", help="end-to-end validation on simulated cohorts")
    s.add_argument("--out", required=True)
    s.add_argument("--panel-size", dest="panel_size", type=int, default=50)
    s.add_argument("--enrich", type=float, default=50.0)
    s.add_argument("--depth", type=int, default=5000)
    s.add_argument("--batch-sd", dest="batch_sd", type=float, default=0.4,
                   help="per-sample shared error multiplier sd (cross-site corr)")
    s.add_argument("--n-train", dest="n_train", type=int, default=200)
    s.add_argument("--n-blank", dest="n_blank", type=int, default=200)
    s.add_argument("--n-rep", dest="n_rep", type=int, default=24)
    s.add_argument("--vaf-levels", dest="vaf_levels", type=float, nargs="+",
                   default=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
    s.add_argument("--alpha", type=float, default=0.05)
    s.add_argument("--null", choices=["empirical", "independent"],
                   default="empirical",
                   help="empirical (covariance-preserving) vs the old independent null")
    s.add_argument("--min-molecules", dest="min_molecules", type=int, default=1000)
    s.add_argument("--iterations", type=int, default=10000)
    s.add_argument("--bootstrap", type=int, default=200)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=cmd_simulate)

    r = sub.add_parser("run-real", help="evaluate real samples from a manifest TSV")
    r.add_argument("--manifest", required=True)
    r.add_argument("--background", required=True)
    r.add_argument("--pon", required=False,
                   help="panel-of-normals TSV for the empirical null")
    r.add_argument("--out", required=True)
    r.add_argument("--alpha", type=float, default=0.05)
    r.add_argument("--min-molecules", dest="min_molecules", type=int, default=1000)
    r.add_argument("--iterations", type=int, default=10000)
    r.add_argument("--bootstrap", type=int, default=200)
    r.add_argument("--seed", type=int, default=0)
    r.set_defaults(func=cmd_run_real)

    t = sub.add_parser("selftest", help="assert mechanics + null-inflation detection")
    t.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
