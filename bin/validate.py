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
# simulate-from-real: a semi-synthetic 2Strands cohort.
#   Real ILM2 blanks give the per-site DEPTH + ERROR structure (the realism over
#   pure simulation); a synthetic forward enrichment + dilution plants the tumor
#   signal. Honest scope: this validates the integration (background -> calibrate
#   -> de-bias -> call) on real coverage/noise structure. It does NOT validate the
#   enrichment MODEL -- that is circular unless the forward model differs from the
#   de-biaser's. The `saturating` knob makes it differ on purpose (robustness).
# ---------------------------------------------------------------------------
@dataclass
class ProbeR:
    site_id: str
    error_rate: float   # real per-site background error (from ILM2 blanks)
    enrich: float       # TRUE assay enrichment factor E (odds multiplier)
    depth: int          # realistic per-site depth (real structure x target_depth)


def enriched_vaf(v, E, model="constant", sat=0.0):
    """Forward enrichment in odds-space -> observed mutant-allele fraction.
      constant  : odds_obs = E * odds_true                    (what the de-biaser assumes)
      saturating: odds_obs = E*odds_true / (1 + sat*odds_true) (effective enrichment
                  wanes as input rises -- agrees with constant as v->0, diverges
                  higher up; sat=0 reduces to constant).
    """
    if v <= 0:
        return 0.0
    o = v / (1.0 - v)
    o_obs = (E * o) / (1.0 + sat * o) if model == "saturating" else E * o
    return o_obs / (1.0 + o_obs)


def load_template(paths):
    """Real blank site-count TSVs -> (depths, error_rates) arrays (interrogate.py
    columns, with plain depth/alt fallbacks)."""
    dep, err = [], []
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                d = int(row.get("depth_unique_molecules") or row.get("depth") or 0)
                a = int(row.get("alt_molecule_count") or row.get("alt") or 0)
                if d > 0:
                    dep.append(d)
                    err.append(min(0.5, a / d))
    return np.array(dep, float), np.array(err, float)


def synth_template(rng, n=400, mean_depth=200, mean_err=5e-4):
    """Fallback template when no real ILM2 counts are supplied (selftest/demo)."""
    dep = rng.lognormal(np.log(mean_depth), 0.5, n)
    err = np.minimum(0.5, rng.lognormal(np.log(mean_err), 0.5, n))
    return dep, err


def make_panel_from_template(rng, n, E, dep, err, target_depth):
    """Bootstrap n panel sites' (error_rate, relative depth) from the template;
    scale the real coverage non-uniformity to an absolute `target_depth`."""
    idx = rng.integers(0, len(err), size=n)
    mean_d = max(1.0, float(np.mean(dep)))
    return [ProbeR(f"site{i}", float(err[j]), E,
                   max(1, int(round(target_depth * dep[j] / mean_d))))
            for i, j in enumerate(idx)]


def simulate_counts_r(rng, panel, vaf, batch_sd, model, sat):
    """One sample's per-site (depth, alt), real per-site depth + forward model."""
    m = float(np.exp(rng.normal(0.0, batch_sd))) if batch_sd > 0 else 1.0
    out = {}
    for p in panel:
        e = int(rng.binomial(p.depth, min(0.5, p.error_rate * m)))
        t = int(rng.binomial(p.depth, enriched_vaf(vaf, p.enrich, model, sat))) if vaf > 0 else 0
        out[p.site_id] = (p.depth, min(p.depth, e + t))
    return out


def build_background_r(rng, panel, n_train, batch_sd):
    alt = {p.site_id: [] for p in panel}
    dep = {p.site_id: [] for p in panel}
    for _ in range(n_train):                       # blanks: vaf=0, model irrelevant
        for sid, (d, a) in simulate_counts_r(rng, panel, 0.0, batch_sd, "constant", 0.0).items():
            alt[sid].append(a)
            dep[sid].append(d)
    bg, donor = {}, {}
    for p in panel:
        av, dv = np.array(alt[p.site_id]), np.array(dep[p.site_id])
        bg[p.site_id] = fit_betabinom_mom(av, dv)
        donor[(p.site_id, 0, "C", "T")] = av / np.maximum(dv, 1)
    return bg, donor


def build_enrichment_r(rng, panel, bg, ladder, reps, model, sat):
    """Calibrate per-probe enrichment from the dilution ladder via the FORWARD
    model. Under `saturating` the engine fits a constant E to non-constant data
    -> a deliberately biased calibration (the robustness test)."""
    enrich = {}
    for p in panel:
        bg_rate = bg[p.site_id][0] / (bg[p.site_id][0] + bg[p.site_id][1])
        known, obs = [], []
        for v in ladder:
            for _ in range(reps):
                _, a = simulate_counts_r(rng, [p], v, 0.0, model, sat)[p.site_id]
                known.append(v)
                obs.append(max(0.0, (a - p.depth * bg_rate) / p.depth))
        enrich[p.site_id] = fit_enrichment_from_dilution(np.array(known), np.array(obs))
    return enrich


def run_from_real(cfg, rng):
    dep, err = cfg["template"]
    model, sat, E = cfg["model"], cfg["sat"], cfg["enrich"]
    panel = make_panel_from_template(rng, cfg["panel_size"], E, dep, err, cfg["target_depth"])
    bg, donor = build_background_r(rng, panel, cfg["n_train"], cfg["batch_sd"])
    enrich = build_enrichment_r(rng, panel, bg, cfg["ladder"], cfg.get("enrich_reps", 5), model, sat)
    fitted_E = float(np.mean(list(enrich.values())))
    dr = donor if cfg.get("use_empirical", True) else None

    blanks = [counts_to_sites(simulate_counts_r(rng, panel, 0.0, cfg["batch_sd"], model, sat), bg, enrich)
              for _ in range(cfg["n_blank"])]
    bc, btf = eval_samples(blanks, rng, cfg["alpha"], cfg["min_molecules"],
                           cfg["n_iter"], cfg["n_boot"], donor_rates=dr)
    fpr, lob = lob_and_fpr(bc, btf)

    level_results = []
    for v in cfg["ladder"]:
        sl = [counts_to_sites(simulate_counts_r(rng, panel, v, cfg["batch_sd"], model, sat), bg, enrich)
              for _ in range(cfg["n_rep"])]
        calls, tfs = eval_samples(sl, rng, cfg["alpha"], cfg["min_molecules"],
                                  cfg["n_iter"], cfg["n_boot"], donor_rates=dr)
        level_results.append((v, calls, tfs))
    table, lod95, loq = lod_loq(level_results)

    pv = cfg["patient_vaf"]                          # the dummy patient (held out)
    psites = counts_to_sites(simulate_counts_r(rng, panel, pv, cfg["batch_sd"], model, sat), bg, enrich)
    pr = call_mrd(psites, rng, p_threshold=cfg["alpha"], min_total_molecules=cfg["min_molecules"],
                  n_iter=cfg["n_iter"], n_boot=cfg["n_boot"], donor_rates=dr)
    rec = pr["tumor_fraction"]
    return {
        "mode": "simulate-from-real", "enrichment_model": model, "saturation": sat,
        "template_sites": int(len(dep)), "target_depth": cfg["target_depth"],
        "enrichment_true": E, "enrichment_fitted": round(fitted_E, 2),
        "enrichment_recovery_ratio": round(fitted_E / E, 3),
        "patient_true_vaf": pv, "patient_recovered_tf": rec, "patient_call": pr["call"],
        "patient_vaf_bias": round(abs(rec - pv) / pv, 3) if pv > 0 else None,
        "nominal_alpha": cfg["alpha"], "observed_fpr": round(fpr, 4),
        "specificity": round(1 - fpr, 4), "limit_of_blank_tf": lob,
        "lod95_vaf": lod95, "loq_vaf": loq, "dilution": table,
        "specificity_holds": bool(fpr <= cfg["alpha"] + 0.02),
        "null_model": "empirical_donor_resampling" if dr is not None else "independent_betabinom",
    }


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


def cmd_simulate_from_real(args):
    rng = np.random.default_rng(args.seed)
    if args.template_counts:
        dep, err = load_template(args.template_counts)
        if len(dep) == 0:
            sys.exit("ERROR: no usable sites in --template-counts")
        src = f"real ({len(dep)} sites from {len(args.template_counts)} blank TSVs)"
    else:
        dep, err = synth_template(rng)
        src = "synthetic (no --template-counts given)"
    cfg = dict(template=(dep, err), panel_size=args.panel_size, enrich=args.enrich,
               model=args.enrichment_model, sat=args.saturation,
               target_depth=args.target_depth, batch_sd=args.batch_sd,
               n_train=args.n_train, n_blank=args.n_blank, n_rep=args.n_rep,
               ladder=args.vaf_ladder, patient_vaf=args.patient_vaf, alpha=args.alpha,
               min_molecules=args.min_molecules, n_iter=args.iterations,
               n_boot=args.bootstrap, use_empirical=args.null == "empirical")
    r = run_from_real(cfg, rng)
    write_report(r, args.out)
    print(f"semi-synthetic 2Strands validation report -> {args.out}")
    print(f"  template          : {src}  target_depth={args.target_depth}")
    print(f"  enrichment model  : {r['enrichment_model']} (saturation={r['saturation']})")
    print(f"  enrichment        : true {r['enrichment_true']}  fitted "
          f"{r['enrichment_fitted']}  recovery {r['enrichment_recovery_ratio']}x")
    print(f"  dummy patient     : true VAF {r['patient_true_vaf']:.2e}  recovered TF "
          f"{r['patient_recovered_tf']:.2e}  bias {r['patient_vaf_bias']}  call {r['patient_call']}")
    print(f"  observed FPR      : {r['observed_fpr']} "
          f"({'OK' if r['specificity_holds'] else 'INFLATED'})  "
          f"LoD95 {r['lod95_vaf']}  LoQ {r['loq_vaf']}")


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

    # --- simulate-from-real: de-biasing on real-template structure, and the
    #     saturating-model robustness knob must actually degrade enrichment recovery ---
    dep, err = synth_template(np.random.default_rng(11), n=150)
    frc = dict(template=(dep, err), panel_size=20, enrich=20.0, target_depth=4000,
               batch_sd=0.2, n_train=60, n_blank=60, n_rep=8, ladder=[1e-2, 1e-3, 1e-4],
               patient_vaf=5e-4, alpha=0.05, min_molecules=1000, n_iter=800, n_boot=30,
               use_empirical=True)
    r_const = run_from_real({**frc, "model": "constant", "sat": 0.0}, np.random.default_rng(3))
    r_sat = run_from_real({**frc, "model": "saturating", "sat": 50.0}, np.random.default_rng(3))
    assert 0.7 <= r_const["enrichment_recovery_ratio"] <= 1.3, r_const  # well-specified -> ~1.0
    assert r_sat["enrichment_recovery_ratio"] < r_const["enrichment_recovery_ratio"], \
        (r_const["enrichment_recovery_ratio"], r_sat["enrichment_recovery_ratio"])  # knob bites

    print("validate.py selftest: PASS")
    print(f"  from-real constant : enrich recovery {r_const['enrichment_recovery_ratio']}x, "
          f"patient bias {r_const['patient_vaf_bias']}")
    print(f"  from-real saturating: enrich recovery {r_sat['enrichment_recovery_ratio']}x "
          f"(misspecified -> recovery degrades, as intended)")
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

    fr = sub.add_parser("simulate-from-real",
                        help="semi-synthetic 2Strands cohort from a real ILM2 depth/error template")
    fr.add_argument("--out", required=True)
    fr.add_argument("--template-counts", dest="template_counts", nargs="*",
                    help="real BLANK site-count TSVs (interrogate.py output) as the "
                         "depth/error template; omit for a synthetic template")
    fr.add_argument("--enrich", type=float, default=20.0,
                    help="TRUE enrichment factor E (odds multiplier)")
    fr.add_argument("--enrichment-model", dest="enrichment_model",
                    choices=["constant", "saturating"], default="constant")
    fr.add_argument("--saturation", type=float, default=0.0,
                    help="saturating-model knob (0=constant; >0 = enrichment wanes with input)")
    fr.add_argument("--target-depth", dest="target_depth", type=int, default=3000,
                    help="absolute per-site depth the real coverage structure is scaled to")
    fr.add_argument("--vaf-ladder", dest="vaf_ladder", type=float, nargs="+",
                    default=[1e-2, 1e-3, 3e-4, 1e-4])
    fr.add_argument("--patient-vaf", dest="patient_vaf", type=float, default=5e-4,
                    help="the dummy patient's true VAF (held out from calibration)")
    fr.add_argument("--panel-size", dest="panel_size", type=int, default=50)
    fr.add_argument("--batch-sd", dest="batch_sd", type=float, default=0.3)
    fr.add_argument("--n-train", dest="n_train", type=int, default=150)
    fr.add_argument("--n-blank", dest="n_blank", type=int, default=120)
    fr.add_argument("--n-rep", dest="n_rep", type=int, default=20)
    fr.add_argument("--null", choices=["empirical", "independent"], default="empirical")
    fr.add_argument("--alpha", type=float, default=0.05)
    fr.add_argument("--min-molecules", dest="min_molecules", type=int, default=1000)
    fr.add_argument("--iterations", type=int, default=10000)
    fr.add_argument("--bootstrap", type=int, default=200)
    fr.add_argument("--seed", type=int, default=0)
    fr.set_defaults(func=cmd_simulate_from_real)

    t = sub.add_parser("selftest", help="assert mechanics + null-inflation detection")
    t.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
