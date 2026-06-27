#!/usr/bin/env python3
"""build_background.py — fit the empirical null from tumor-free (blank) runs.

Turns blank cfDNA interrogate.py outputs into the two files mrd_integrate needs:

  background.tsv : chrom pos ref alt alpha beta enrich    per-site beta-binomial
  pon.tsv        : chrom pos ref alt <donor1> <donor2> ...  per-donor error rates
                   (the covariance-preserving empirical null)

Inputs are blank site-count TSVs (interrogate.py output: columns
depth_unique_molecules, alt_molecule_count, ...). KEEP FOLDS DISJOINT: build the
background from blanks that are NOT in your run-real specificity test set
(see VALIDATION.md Step 3).

Run:
  build_background.py run --out-background background.tsv --out-pon pon.tsv \
      results/interrogate/SRR1.tsv results/interrogate/SRR2.tsv ...
  # or pull the blank rows straight from a run-real manifest:
  build_background.py run --manifest train_blanks.manifest.tsv ...
Self-test:
  build_background.py selftest
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mrd_integrate import (fit_betabinom_mom, load_sites,                # noqa: E402
                           load_donor_rates, call_mrd)

Key = tuple  # (chrom, pos:int, ref, alt)


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------
def read_counts(path: str) -> dict:
    """One blank's interrogate.py output -> {site key: (depth, alt)}."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
            out[k] = (int(r["depth_unique_molecules"]), int(r["alt_molecule_count"]))
    if not out:
        sys.exit(f"ERROR: no site rows in {path}")
    return out


def load_enrich(path: str | None) -> dict:
    """Optional per-probe enrichment (chrom,pos,ref,alt,enrich). Default 1.0."""
    if not path:
        return {}
    return {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]): float(r["enrich"])
            for r in csv.DictReader(open(path, newline=""), delimiter="\t")}


def blanks_from_manifest(manifest: str, label: str) -> list[str]:
    paths = []
    with open(manifest, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("label") == label:
                paths.append(r["site_counts_path"])
    return paths


# ---------------------------------------------------------------------------
# Core (pure, testable)
# ---------------------------------------------------------------------------
def build(paths: list[str], enrich_map: dict | None = None):
    """Returns (bg_rows, pon_rows, donor_names).

    Sites = the intersection across all blanks (interrogate emits every panel
    site per run, so this normally equals the full panel; we warn if not).
    Per site: beta-binomial (alpha,beta) by method-of-moments over donors with
    depth>0, and the per-donor error-rate vector for the empirical null.
    """
    if not paths:
        sys.exit("ERROR: no blank count files given.")
    per = [read_counts(p) for p in paths]
    sets = [set(d) for d in per]
    common, allk = set.intersection(*sets), set.union(*sets)
    if common != allk:
        print(f"WARN: {len(allk - common)} site(s) absent from some of the "
              f"{len(paths)} blanks; using the intersection ({len(common)} sites).",
              file=sys.stderr)
    sites = sorted(common, key=lambda k: (k[0], k[1], k[2], k[3]))
    if not sites:
        sys.exit("ERROR: blanks share no common sites — are these the same panel?")

    enrich_map = enrich_map or {}
    bg_rows, pon_rows = [], []
    for k in sites:
        depths = np.array([per[i][k][0] for i in range(len(paths))])
        alts = np.array([per[i][k][1] for i in range(len(paths))])
        keep = depths > 0
        if keep.sum() < 1:
            print(f"WARN: site {k} has zero depth in all blanks; skipped.",
                  file=sys.stderr)
            continue
        alpha, beta = fit_betabinom_mom(alts[keep], depths[keep])
        bg_rows.append((*k, alpha, beta, enrich_map.get(k, 1.0)))
        rates = np.where(depths > 0, alts / np.maximum(depths, 1), 0.0)
        pon_rows.append((*k, rates))
    return bg_rows, pon_rows, [_donor_name(p) for p in paths]


def _donor_name(path: str) -> str:
    b = os.path.basename(path)
    return b[:-4] if b.endswith(".tsv") else b


def write_background(rows, path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "alpha", "beta", "enrich"])
        for c, p, r, a, al, be, en in rows:
            w.writerow([c, p, r, a, f"{al:.8g}", f"{be:.8g}", f"{en:.8g}"])


def write_pon(rows, donor_names, path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", *donor_names])
        for c, p, r, a, rates in rows:
            w.writerow([c, p, r, a, *[f"{x:.8g}" for x in rates]])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def cmd_run(args):
    paths = list(args.counts)
    if args.manifest:
        paths += blanks_from_manifest(args.manifest, args.label)
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        sys.exit(f"ERROR: {len(missing)} count file(s) not found, e.g. {missing[0]}")
    bg, pon, donors = build(paths, load_enrich(args.enrich_file))
    write_background(bg, args.out_background)
    write_pon(pon, donors, args.out_pon)
    print(f"background: {len(bg)} sites from {len(donors)} blank donors "
          f"-> {args.out_background}, {args.out_pon}")
    if len(donors) < 4:
        print(f"WARN: only {len(donors)} blank donors; the empirical null is "
              f"granular. More blanks = a better-resolved null tail.", file=sys.stderr)


def cmd_selftest(_args):
    import tempfile, shutil
    rng = np.random.default_rng(0)
    d = tempfile.mkdtemp()
    cols = ["chrom", "pos", "ref", "alt", "depth_unique_molecules",
            "alt_molecule_count", "duplex_support", "mean_consensus_qual"]

    def write_counts(path, sites, depth, err):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
            for (c, p, r, a) in sites:
                w.writerow([c, p, r, a, depth, int(rng.binomial(depth, err)), 0, 35.0])

    try:
        sites = [("chr1", 1000 + i, "C", "T") for i in range(8)]
        depth, err, ndonor = 5000, 5e-4, 6
        blanks = [os.path.join(d, f"blank{j}.tsv") for j in range(ndonor)]
        for p in blanks:
            write_counts(p, sites, depth, err)

        bg, pon, donors = build(blanks)
        bgp, ponp = os.path.join(d, "background.tsv"), os.path.join(d, "pon.tsv")
        write_background(bg, bgp); write_pon(pon, donors, ponp)

        # 1) fitted error rate recovers the truth
        mean_err = float(np.mean([al / (al + be) for (_, _, _, _, al, be, _) in bg]))
        assert 1e-4 < mean_err < 2e-3, mean_err
        # 2) the files load back into the engine's readers with the right shape
        dr = load_donor_rates(ponp)
        assert len(dr) == len(sites), len(dr)
        assert all(len(v) == ndonor for v in dr.values())
        # 3) round-trip: a fresh blank sample scores negative through the engine,
        #    using the background + empirical null we just built
        smp_path = os.path.join(d, "sample.tsv")
        write_counts(smp_path, sites, depth, err)
        sites_obj = load_sites(smp_path, bgp)
        res = call_mrd(sites_obj, np.random.default_rng(1), n_iter=2000, n_boot=100,
                       min_total_molecules=100, donor_rates=dr)
        assert res["call"] in ("negative", "indeterminate"), res
        assert res["null_model"] == "empirical_donor_resampling", res

        print("build_background.py selftest: PASS "
              f"(sites={len(sites)} donors={ndonor} mean_err={mean_err:.2e} "
              f"blank_call={res['call']} p={res['pvalue']:.3f})")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="build background.tsv + pon.tsv from blanks")
    r.add_argument("counts", nargs="*", help="blank site-count TSVs (interrogate output)")
    r.add_argument("--manifest", help="run-real manifest; uses rows with --label")
    r.add_argument("--label", default="blank", help="manifest label to treat as blank")
    r.add_argument("--enrich-file", dest="enrich_file",
                   help="optional per-probe enrichment TSV (chrom,pos,ref,alt,enrich)")
    r.add_argument("--out-background", dest="out_background", default="background.tsv")
    r.add_argument("--out-pon", dest="out_pon", default="pon.tsv")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("selftest", help="round-trip through the engine on synthetic blanks")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
