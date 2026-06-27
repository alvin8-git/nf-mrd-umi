#!/usr/bin/env python3
"""mrd_integrate.py — Pipeline B / Stage B6: the MRD LoD engine (Approach A).

Tumor-informed, panel-integrated MRD detection. Given per-site molecule counts
(from interrogate.py) and a per-site beta-binomial error null (from a frozen
healthy-donor cfDNA cohort), compute a single patient-level call:

    tumor_fraction + Monte-Carlo p-value + bootstrap CI + QC/LoD floor.

Detection is a COMPOSITE call across all panel sites, never per-variant — that
is how sub-0.1% is reached when no single locus has enough molecules (design P2).

Null per site i:  alt molecules a_i ~ BetaBinomial(n=d_i, alpha_i, beta_i),
where (alpha_i, beta_i) describe that locus's background error-rate distribution
across healthy donors. Beta-binomial (not plain binomial) absorbs cross-
donor/run overdispersion, which is the realistic failure mode.

Run:       mrd_integrate.py run --site-counts S.tsv --background B.tsv \
                                --panel panel.vcf --out sample.mrd.json
Self-test: mrd_integrate.py selftest
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass

import numpy as np
from scipy.stats import betabinom


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
@dataclass
class Site:
    chrom: str
    pos: int
    ref: str
    alt: str
    depth: int          # d_i : unique molecules covering the site
    alt_count: int      # a_i : molecules supporting the tumor's ALT allele
    alpha: float        # beta-binomial null shape a (healthy-donor cohort)
    beta: float         # beta-binomial null shape b
    enrich: float = 1.0 # per-probe allele-specific enrichment factor: how much
                        # the 2Strands probe over-represents the MUTANT allele vs
                        # unbiased sampling. 1.0 = none. Calibrated from a
                        # contrived dilution series. Used ONLY to de-bias the
                        # reported tumor fraction — NOT the detection p-value,
                        # which stays in enriched space vs an enriched null.

    @property
    def error_rate(self) -> float:
        return self.alpha / (self.alpha + self.beta)


# ---------------------------------------------------------------------------
# Core statistics (pure, unit-testable)
# ---------------------------------------------------------------------------
def estimate_vaf(sites: list[Site], apply_enrichment: bool = True) -> float:
    """Background-subtracted tumor fraction across the panel.

    Per site: observed excess VAF = max(0, (a_i - d_i*e_i)/d_i). Under 2Strands
    allele-specific enrichment the probe over-represents the mutant allele by a
    per-probe factor enrich_i, so the OBSERVED excess is inflated; we de-bias by
    dividing each site's excess by enrich_i (low-VAF linear approximation; the
    exact odds-space correction is calibrated from the dilution series). Sites
    are depth-weighted. apply_enrichment=False returns the raw observed
    (enriched) VAF for QC.

    Reported as tumor fraction (VAF). For het SNVs cellular fraction ~= 2*VAF;
    het conversion is left to the clinical layer (zygosity/CN aware).
    """
    total_depth = sum(s.depth for s in sites)
    if total_depth == 0:
        return 0.0
    weighted = 0.0
    for s in sites:
        if s.depth == 0:
            continue
        excess_frac = max(0.0, (s.alt_count - s.depth * s.error_rate) / s.depth)
        e = s.enrich if (apply_enrichment and s.enrich > 0) else 1.0
        # De-enrich in ODDS space (enrichment multiplies the mutant-allele odds):
        # o_true = o_obs / e ; vaf_true = o/(1+o). At low fractions this ~= /e,
        # but it stays accurate at the high OBSERVED fractions enrichment creates.
        f = min(excess_frac, 1.0 - 1e-9)
        o_true = (f / (1.0 - f)) / e
        weighted += s.depth * (o_true / (1.0 + o_true))
    return weighted / total_depth


def fit_enrichment_from_dilution(known_vaf: np.ndarray,
                                 observed_vaf: np.ndarray) -> float:
    """Per-probe enrichment factor from a contrived dilution series.

    Allele-specific enrichment multiplies the mutant-allele ODDS:
        odds_obs = enrich * odds_true,   odds = vaf / (1 - vaf).
    So enrich is the through-origin slope in ODDS space. Unlike a slope in VAF
    space, this does NOT bias when high calibration points saturate toward
    VAF=1 (the old bug: fit recovered ~0.43x of truth -> TF over-estimated
    ~2.3x). Points at 0/1 are dropped; result floored at 1.0.
    """
    known = np.asarray(known_vaf, float)
    obs = np.asarray(observed_vaf, float)
    ok = (known > 0) & (known < 1) & (obs > 0) & (obs < 1)
    if not np.any(ok):
        return 1.0
    o_known = known[ok] / (1.0 - known[ok])
    o_obs = obs[ok] / (1.0 - obs[ok])
    denom = float(np.sum(o_known * o_known))
    if denom <= 0:
        return 1.0
    return max(1.0, float(np.sum(o_known * o_obs) / denom))


def monte_carlo_pvalue(sites: list[Site], rng: np.random.Generator,
                       n_iter: int = 20000, donor_rates: dict | None = None
                       ) -> float:
    """One-sided MC p-value for total alt signal S = sum_i a_i, depth-matched.

    Two nulls:
    - EMPIRICAL (donor_rates given) — the correct one. For each replicate, draw a
      WHOLE healthy donor's per-site error-rate vector and sample alt_i ~
      Binomial(d_i, rate_i^donor). Because a single donor's rates are taken as a
      unit, cross-site correlation (batch/deamination/GC moving all sites
      together) is preserved in the null tail. This is what stops a globally-hot
      sample from reading as positive. donor_rates maps site key -> 1D array of
      per-donor error rates.
    - INDEPENDENT (donor_rates None) — fallback only. Per-site beta-binomial
      draws, summed assuming independence. ANTICONSERVATIVE under correlated
      error (validation harness measures ~5.7x inflated FPR). Kept for back-compat.

    +1 smoothing keeps p bounded away from 0 (never claim p=0 from finite draws).
    """
    informative = [s for s in sites if s.depth > 0]
    if not informative:
        return 1.0
    s_obs = sum(s.alt_count for s in informative)

    if donor_rates is not None:
        keys = [(s.chrom, s.pos, s.ref, s.alt) for s in informative]
        missing = [k for k in keys if k not in donor_rates]
        if missing:
            sys.exit(f"ERROR: empirical null missing donor rates for {missing[:3]}"
                     f"{'...' if len(missing) > 3 else ''}")
        R = np.array([donor_rates[k] for k in keys])      # [n_sites, n_donors]
        d = np.array([s.depth for s in informative])      # [n_sites]
        donor_idx = rng.integers(0, R.shape[1], size=n_iter)
        rates = R[:, donor_idx].T                         # [n_iter, n_sites]
        null_totals = rng.binomial(d, rates).sum(axis=1)  # d broadcasts per row
    else:
        null_totals = np.zeros(n_iter, dtype=np.int64)
        for s in informative:
            null_totals += betabinom.rvs(s.depth, s.alpha, s.beta,
                                         size=n_iter, random_state=rng)

    n_ge = int(np.count_nonzero(null_totals >= s_obs))
    return (1 + n_ge) / (n_iter + 1)


def bootstrap_ci(sites: list[Site], rng: np.random.Generator,
                 n_boot: int = 2000, alpha: float = 0.05,
                 apply_enrichment: bool = True) -> tuple[float, float]:
    """Percentile CI for tumor fraction by resampling sites with replacement."""
    informative = [s for s in sites if s.depth > 0]
    if len(informative) < 2:
        v = estimate_vaf(informative, apply_enrichment)
        return (v, v)
    idx = np.arange(len(informative))
    vafs = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(idx, size=len(idx), replace=True)
        vafs[b] = estimate_vaf([informative[i] for i in pick], apply_enrichment)
    lo = float(np.percentile(vafs, 100 * alpha / 2))
    hi = float(np.percentile(vafs, 100 * (1 - alpha / 2)))
    return (lo, hi)


def fit_betabinom_mom(donor_alt: np.ndarray, donor_depth: np.ndarray,
                      min_kappa: float = 2.0) -> tuple[float, float]:
    """Method-of-moments beta-binomial fit for ONE site from a healthy cohort.

    donor_alt/donor_depth: per-donor alt and total molecule counts at this site.
    Returns (alpha, beta). Falls back to a tight pseudo-count prior when the
    cohort shows no variance (common for clean sites), so the null is never
    degenerate. Used offline to build the background file; exposed here so the
    same code path is what gets unit-tested.
    """
    p = donor_alt / np.maximum(donor_depth, 1)
    m = float(np.mean(p))
    v = float(np.var(p))
    if m <= 0:
        # No errors observed: tight prior at the cohort's resolution limit.
        n_eff = float(np.sum(donor_depth)) or 1.0
        m = 1.0 / (n_eff + 1.0)
        v = m * (1 - m) / (n_eff + 1.0)
    if v <= 0 or m * (1 - m) <= v:
        kappa = max(min_kappa, np.sum(donor_depth) / max(len(donor_depth), 1))
    else:
        kappa = max(min_kappa, m * (1 - m) / v - 1.0)
    alpha = max(m * kappa, 1e-6)
    beta = max((1 - m) * kappa, 1e-6)
    return (alpha, beta)


def call_mrd(sites: list[Site], rng: np.random.Generator,
             p_threshold: float = 0.05,
             min_total_molecules: int = 1000,
             min_informative_sites: int = 1,
             n_iter: int = 20000, n_boot: int = 2000,
             donor_rates: dict | None = None) -> dict:
    """Full panel-integrated MRD call with QC/LoD gating.

    donor_rates (site key -> per-donor error-rate array) selects the
    covariance-preserving EMPIRICAL null; omit it only for the back-compat
    independent null (anticonservative under correlated error)."""
    informative = [s for s in sites if s.depth > 0]
    total_molecules = sum(s.depth for s in informative)
    n_informative = len(informative)

    true_tf = estimate_vaf(informative, apply_enrichment=True)
    observed_vaf = estimate_vaf(informative, apply_enrichment=False)
    pvalue = monte_carlo_pvalue(informative, rng, n_iter=n_iter,
                                donor_rates=donor_rates)
    ci_lo, ci_hi = bootstrap_ci(informative, rng, n_boot=n_boot,
                                apply_enrichment=True)
    enriched = any(s.enrich != 1.0 for s in informative)

    # LoD gate first: too few molecules => we cannot rule tumor in OR out.
    if total_molecules < min_total_molecules or n_informative < min_informative_sites:
        call = "indeterminate"
    elif pvalue < p_threshold:
        call = "positive"
    else:
        call = "negative"

    return {
        "call": call,
        "tumor_fraction": true_tf,                 # enrichment-de-biased plasma VAF
        "observed_vaf_enriched": observed_vaf,     # pre-correction, for QC
        "enrichment_corrected": enriched,
        "tumor_fraction_ci95": [ci_lo, ci_hi],
        "pvalue": pvalue,
        "pvalue_space": "enriched (observed counts vs enriched null)",
        "null_model": "empirical_donor_resampling" if donor_rates is not None
                      else "independent_betabinom",
        "p_threshold": p_threshold,
        "n_sites_total": len(sites),
        "n_sites_informative": n_informative,
        "unique_molecular_depth": total_molecules,
        "lod_min_molecules": min_total_molecules,
        "alt_molecules_observed": sum(s.alt_count for s in informative),
        "alt_molecules_expected_bg": round(
            sum(s.depth * s.error_rate for s in informative), 3),
    }


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def _read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_sites(site_counts_tsv: str, background_tsv: str) -> list[Site]:
    """Join interrogate.py counts with the per-site beta-binomial background.

    Both keyed on (chrom,pos,ref,alt). A panel site with no background entry is
    a hard error: an uncharacterized locus cannot be scored (fail loud).
    """
    # background carries (alpha, beta, enrich); enrich optional, defaults to 1.0
    bg = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]):
          (float(r["alpha"]), float(r["beta"]), float(r.get("enrich", 1.0) or 1.0))
          for r in _read_tsv(background_tsv)}
    sites: list[Site] = []
    for r in _read_tsv(site_counts_tsv):
        key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
        if key not in bg:
            sys.exit(f"ERROR: no background for panel site {key}; "
                     f"rebuild background or drop the site.")
        alpha, beta, enrich = bg[key]
        d = int(r["depth_unique_molecules"])
        a = int(r["alt_molecule_count"])
        if a > d:
            sys.exit(f"ERROR: alt({a}) > depth({d}) at {key} — bad input.")
        sites.append(Site(r["chrom"], int(r["pos"]), r["ref"], r["alt"],
                          d, a, alpha, beta, enrich))
    if not sites:
        sys.exit("ERROR: no panel sites after join.")
    return sites


def load_donor_rates(pon_tsv: str) -> dict:
    """Healthy-donor panel-of-normals for the empirical null. TSV: chrom,pos,ref,
    alt, then one column per donor holding that donor's per-site error RATE
    (alt/depth). Returns {site key -> 1D array of per-donor rates}."""
    rows = _read_tsv(pon_tsv)
    if not rows:
        sys.exit("ERROR: empty panel-of-normals.")
    meta = {"chrom", "pos", "ref", "alt"}
    donor_cols = [c for c in rows[0].keys() if c not in meta]
    return {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]):
            np.array([float(r[c]) for c in donor_cols]) for r in rows}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def cmd_run(args) -> None:
    rng = np.random.default_rng(args.seed)
    sites = load_sites(args.site_counts, args.background)
    donor_rates = load_donor_rates(args.pon) if args.pon else None
    if donor_rates is None:
        print("WARNING: no --pon panel-of-normals; using the independent null "
              "(anticonservative under correlated error).", file=sys.stderr)
    result = call_mrd(
        sites, rng,
        p_threshold=args.p_threshold,
        min_total_molecules=args.min_molecules,
        min_informative_sites=args.min_sites,
        n_iter=args.iterations, n_boot=args.bootstrap,
        donor_rates=donor_rates,
    )
    result["patient_id"] = args.patient_id
    result["timepoint"] = args.timepoint
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"{result['call']}  TF={result['tumor_fraction']:.2e}  "
          f"(obs {result['observed_vaf_enriched']:.2e})  "
          f"p={result['pvalue']:.2e}  -> {args.out}")


def _make_synthetic(rng, n_sites=50, depth=5000, err=5e-4, inject_vaf=0.0,
                    enrich=1.0):
    """Build a synthetic panel; inject_vaf>0 simulates a tumor signal at true
    plasma VAF inject_vaf, observed after allele-specific enrichment `enrich`."""
    kappa = 5000.0  # low overdispersion
    alpha, beta = err * kappa, (1 - err) * kappa
    sites = []
    for i in range(n_sites):
        bg = int(betabinom.rvs(depth, alpha, beta, random_state=rng))
        if inject_vaf > 0:                       # odds-model enrichment
            o = inject_vaf / (1.0 - inject_vaf)
            obs_p = (enrich * o) / (1.0 + enrich * o)
            tumor = int(rng.binomial(depth, obs_p))
        else:
            tumor = 0
        a = min(depth, bg + tumor)
        sites.append(Site("chr1", 1000 + i, "C", "T", depth, a, alpha, beta,
                          enrich))
    return sites


def cmd_selftest(_args) -> None:
    rng = np.random.default_rng(12345)

    # MoM fit recovers a sane error rate.
    donor_depth = np.full(40, 5000)
    donor_alt = rng.binomial(5000, 5e-4, size=40)
    a, b = fit_betabinom_mom(donor_alt, donor_depth)
    assert 1e-4 < a / (a + b) < 2e-3, f"MoM error rate off: {a/(a+b)}"

    # Negative sample: pure background -> not significant, TF ~ 0.
    neg = _make_synthetic(rng, inject_vaf=0.0)
    rneg = call_mrd(neg, rng, n_iter=5000, n_boot=500)
    assert rneg["call"] == "negative", rneg
    assert rneg["pvalue"] > 0.05, rneg
    assert rneg["tumor_fraction"] < 5e-4, rneg

    # Positive sample (no enrichment): inject 0.2% -> significant, estimate ~truth.
    pos = _make_synthetic(rng, inject_vaf=2e-3)
    rpos = call_mrd(pos, rng, n_iter=5000, n_boot=500)
    assert rpos["call"] == "positive", rpos
    assert rpos["pvalue"] < 0.01, rpos
    assert abs(rpos["tumor_fraction"] - 2e-3) < 1e-3, rpos

    # Allele-specific enrichment: true 0.2% observed at ~50x -> ~10% enriched,
    # but de-biased tumor_fraction must recover ~0.2%.
    enr = _make_synthetic(rng, inject_vaf=2e-3, enrich=50.0)
    renr = call_mrd(enr, rng, n_iter=5000, n_boot=500)
    assert renr["call"] == "positive", renr
    assert renr["enrichment_corrected"] is True, renr
    assert renr["observed_vaf_enriched"] > 0.05, renr           # inflated
    assert abs(renr["tumor_fraction"] - 2e-3) < 1e-3, renr      # recovered

    # Enrichment factor recovered from a dilution series (odds model), including
    # a high point that saturates VAF -> would have broken the old VAF-space fit.
    known = np.array([1e-3, 5e-3, 1e-2, 5e-2])
    e_true = 47.0
    ok = known / (1.0 - known)
    obs = (e_true * ok) / (1.0 + e_true * ok)
    e_hat = fit_enrichment_from_dilution(known, obs)
    assert abs(e_hat - 47.0) < 1.0, e_hat

    # LoD gate: starve molecular depth -> indeterminate, not negative.
    thin = _make_synthetic(rng, n_sites=2, depth=50, inject_vaf=0.0)
    rthin = call_mrd(thin, rng, min_total_molecules=1000, n_iter=2000, n_boot=200)
    assert rthin["call"] == "indeterminate", rthin

    print("mrd_integrate.py selftest: PASS "
          f"(neg p={rneg['pvalue']:.3f}, pos TF={rpos['tumor_fraction']:.2e}, "
          f"enriched obs={renr['observed_vaf_enriched']:.2e} -> "
          f"TF={renr['tumor_fraction']:.2e}, e_hat={e_hat:.1f})")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="score one cfDNA timepoint")
    r.add_argument("--site-counts", dest="site_counts", required=True)
    r.add_argument("--background", required=True)
    r.add_argument("--pon", required=False,
                   help="healthy-donor panel-of-normals TSV for the empirical "
                        "null (per-site per-donor error rates); strongly advised")
    r.add_argument("--panel", required=False, help="panel VCF (provenance only)")
    r.add_argument("--out", required=True)
    r.add_argument("--patient-id", dest="patient_id", default="NA")
    r.add_argument("--timepoint", default="NA")
    r.add_argument("--p-threshold", dest="p_threshold", type=float, default=0.05)
    r.add_argument("--min-molecules", dest="min_molecules", type=int, default=1000)
    r.add_argument("--min-sites", dest="min_sites", type=int, default=1)
    r.add_argument("--iterations", type=int, default=20000)
    r.add_argument("--bootstrap", type=int, default=2000)
    r.add_argument("--seed", type=int, default=0)
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("selftest", help="run synthetic-data asserts")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
