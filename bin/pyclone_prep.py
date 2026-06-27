#!/usr/bin/env python3
"""pyclone_prep.py — Pipeline A glue around PyClone-vi.

build  : somatic VCF + FACETS copy-number VCF + tumour purity -> PyClone-vi
         input table (one row per somatic SNV).
to-ccf : PyClone-vi results -> the per-variant (chrom pos ref alt ccf clonal_prob)
         table that panel_select.py consumes. The clonal/truncal cluster is the
         one with the highest mean cellular prevalence; clonal_prob is the
         assignment probability for variants in that cluster, else 0.

Run:
  pyclone_prep.py build  --vcf somatic.vcf.gz --facets facets.vcf.gz \
                         --purity purity.txt --out pyclone_input.tsv
  pyclone_prep.py to-ccf --pyclone pyclone.tsv --out ccf.tsv
Self-test:
  pyclone_prep.py selftest
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Pure cores (unit-testable, no pysam)
# ---------------------------------------------------------------------------
def cn_for_site(segments, chrom, pos):
    """segments: list of (chrom, start, end, major_cn, minor_cn). Return the
    (major, minor) of the segment containing (chrom, pos), else diploid (1,1)."""
    for c, s, e, major, minor in segments:
        if c == chrom and s <= pos <= e:
            return major, minor
    return 1, 1   # diploid heterozygous default (major 1 + minor 1 = total 2)


def pyclone_to_ccf(rows):
    """rows: dicts with mutation_id ('chrom:pos:ref:alt'), cluster_id,
    cellular_prevalence, cluster_assignment_prob. Returns dicts
    chrom,pos,ref,alt,ccf,clonal_prob. Truncal cluster = highest mean CCF."""
    sums, counts = defaultdict(float), defaultdict(int)
    for r in rows:
        sums[r["cluster_id"]] += float(r["cellular_prevalence"])
        counts[r["cluster_id"]] += 1
    if not counts:
        return []
    means = {k: sums[k] / counts[k] for k in counts}
    clonal = max(means, key=means.get)

    out = []
    for r in rows:
        parts = r["mutation_id"].split(":")
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts
        ccf = float(r["cellular_prevalence"])
        cp = float(r.get("cluster_assignment_prob", 1.0)) if r["cluster_id"] == clonal else 0.0
        out.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "ccf": round(ccf, 4), "clonal_prob": round(cp, 4)})
    return out


def parse_purity(text: str) -> float:
    """Pull a purity value out of a FACETS purity.txt / VCF header blob."""
    m = re.search(r"purity[=\s:]+([0-9.]+)", text, re.I)
    return float(m.group(1)) if m else 0.5   # fallback; should be reviewed


# ---------------------------------------------------------------------------
# pysam I/O (not exercised by selftest)
# ---------------------------------------------------------------------------
def load_facets_segments(facets_vcf):
    import pysam
    segs = []
    with pysam.VariantFile(facets_vcf) as vf:
        for rec in vf:
            tcn = rec.info.get("TCN_EM")
            if tcn is None:
                continue
            lcn = rec.info.get("LCN_EM") or 0
            segs.append((rec.chrom, rec.pos, int(rec.stop),
                         max(int(tcn) - int(lcn), 0), int(lcn)))
    return segs


def cmd_build(args):
    import pysam
    segs = load_facets_segments(args.facets)
    purity = parse_purity(open(args.purity).read())
    with pysam.VariantFile(args.vcf) as vf, open(args.out, "w", newline="") as out:
        samples = list(vf.header.samples)
        tumor = args.tumor_sample or samples[0]
        if tumor not in samples:
            sys.exit(f"ERROR: tumour sample '{tumor}' not in VCF samples {samples}")
        w = csv.writer(out, delimiter="\t")
        w.writerow(["mutation_id", "sample_id", "ref_counts", "alt_counts",
                    "major_cn", "minor_cn", "normal_cn", "tumour_content"])
        n = 0
        for rec in vf:
            filt = set(rec.filter.keys())
            if filt and "PASS" not in filt:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:   # SNV only
                continue
            ad = rec.samples[tumor].get("AD")
            if not ad or len(ad) < 2:
                continue
            major, minor = cn_for_site(segs, rec.chrom, rec.pos)
            normal_cn = 1 if rec.chrom in ("chrX", "chrY", "X", "Y") else 2
            mid = f"{rec.chrom}:{rec.pos}:{rec.ref}:{rec.alts[0]}"
            w.writerow([mid, tumor, int(ad[0]), int(ad[1]),
                        major, minor, normal_cn, round(purity, 4)])
            n += 1
    print(f"pyclone_prep build: {n} SNVs -> {args.out} (purity {purity:.3f})")


def cmd_to_ccf(args):
    rows = list(csv.DictReader(open(args.pyclone), delimiter="\t"))
    ccf = pyclone_to_ccf(rows)
    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["chrom", "pos", "ref", "alt",
                                            "ccf", "clonal_prob"], delimiter="\t")
        w.writeheader()
        w.writerows(ccf)
    print(f"pyclone_prep to-ccf: {len(ccf)} variants -> {args.out}")


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
def cmd_selftest(_args):
    segs = [("chr1", 100, 200, 2, 1), ("chr2", 50, 80, 1, 0)]
    assert cn_for_site(segs, "chr1", 150) == (2, 1)
    assert cn_for_site(segs, "chr2", 60) == (1, 0)
    assert cn_for_site(segs, "chr3", 10) == (1, 1)     # diploid default
    assert parse_purity("##purity=0.62\n##ploidy=2.0") == 0.62

    rows = [
        {"mutation_id": "chr1:100:C:T", "cluster_id": "0",
         "cellular_prevalence": "0.92", "cluster_assignment_prob": "0.99"},
        {"mutation_id": "chr2:200:G:A", "cluster_id": "0",
         "cellular_prevalence": "0.88", "cluster_assignment_prob": "0.95"},
        {"mutation_id": "chr3:300:A:G", "cluster_id": "1",
         "cellular_prevalence": "0.31", "cluster_assignment_prob": "0.97"},
    ]
    ccf = {r["chrom"]: r for r in pyclone_to_ccf(rows)}
    assert ccf["chr1"]["clonal_prob"] == 0.99, ccf   # truncal cluster
    assert ccf["chr1"]["ccf"] == 0.92, ccf
    assert ccf["chr3"]["clonal_prob"] == 0.0, ccf    # subclonal -> 0
    print("pyclone_prep.py selftest: PASS (cn lookup, purity parse, truncal cluster)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="VCF + FACETS + purity -> PyClone-vi input")
    b.add_argument("--vcf", required=True)
    b.add_argument("--facets", required=True)
    b.add_argument("--purity", required=True)
    b.add_argument("--tumor-sample", dest="tumor_sample")
    b.add_argument("--out", required=True)
    b.set_defaults(func=cmd_build)

    t = sub.add_parser("to-ccf", help="PyClone-vi results -> panel_select ccf TSV")
    t.add_argument("--pyclone", required=True)
    t.add_argument("--out", required=True)
    t.set_defaults(func=cmd_to_ccf)

    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
