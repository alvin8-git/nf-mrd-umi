#!/usr/bin/env python3
"""normal_evidence.py — Pipeline A glue: per-site buffy-coat (normal) support.

For each somatic SNV, pile up the matched buffy-coat (WBC) BAM and report
(normal_alt, normal_depth). panel_select.py uses this to subtract CHIP and
germline (design P4 — CHIP is defined by presence in the WBC compartment, not
by gene). Output columns: chrom pos ref alt normal_alt normal_depth.

Run:       normal_evidence.py pileup --vcf somatic.vcf.gz --bam buffy.bam \
                                      --ref ref.fa --out normal_evidence.tsv
Self-test: normal_evidence.py selftest
"""
from __future__ import annotations

import argparse
import csv

_ACGT = frozenset("ACGT")


# ---------------------------------------------------------------------------
# Pure core (unit-testable, no pysam)
# ---------------------------------------------------------------------------
def count_alt_depth(bases, quals, alt, min_bq: int = 20):
    """Count alt-supporting reads and qualifying depth at one pileup column."""
    if len(bases) != len(quals):
        raise ValueError("bases/quals length mismatch")
    alt = alt.upper()
    depth = alt_ct = 0
    for b, q in zip(bases, quals):
        bb = b.upper()
        if bb not in _ACGT or q < min_bq:
            continue
        depth += 1
        if bb == alt:
            alt_ct += 1
    return alt_ct, depth


# ---------------------------------------------------------------------------
# pysam wrapper (not exercised by selftest)
# ---------------------------------------------------------------------------
def pileup_sites(vcf, bam, out, min_bq: int = 20, min_mapq: int = 20):
    import pysam

    sites = []
    with pysam.VariantFile(vcf) as vf:
        for rec in vf:
            for a in rec.alts or []:
                if len(rec.ref) == 1 and len(a) == 1:   # SNV only
                    sites.append((rec.chrom, rec.pos, rec.ref, a))

    bf = pysam.AlignmentFile(bam, "rb")
    rows = []
    for chrom, pos1, ref, alt in sites:
        pos0 = pos1 - 1
        bases, quals = [], []
        for col in bf.pileup(chrom, pos0, pos0 + 1, truncate=True,
                             min_base_quality=0, stepper="samtools"):
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                if pr.alignment.mapping_quality < min_mapq:
                    continue
                bases.append(pr.alignment.query_sequence[pr.query_position])
                quals.append(pr.alignment.query_qualities[pr.query_position])
        alt_ct, depth = count_alt_depth(bases, quals, alt, min_bq)
        rows.append((chrom, pos1, ref, alt, alt_ct, depth))
    bf.close()

    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "normal_alt", "normal_depth"])
        w.writerows(rows)
    print(f"normal_evidence: {len(rows)} sites -> {out}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def cmd_pileup(args):
    pileup_sites(args.vcf, args.bam, args.out, args.min_bq, args.min_mapq)


def cmd_selftest(_args):
    # 3 ref(C), 1 high-q alt(T), 1 low-q alt(dropped), 1 N(dropped), 1 alt(T)
    bases = ["C", "C", "C", "T", "T", "N", "T"]
    quals = [40,  40,  40,  40,  10,  40,  40]
    a, d = count_alt_depth(bases, quals, "T", min_bq=20)
    assert d == 5 and a == 2, (a, d)              # 5 qualifying, 2 alt

    a0, d0 = count_alt_depth(["C"] * 4, [40] * 4, "T")
    assert a0 == 0 and d0 == 4, (a0, d0)          # no alt support

    try:
        count_alt_depth(["C"], [40, 40], "T")
        raise AssertionError("expected ValueError on length mismatch")
    except ValueError:
        pass
    print("normal_evidence.py selftest: PASS")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("pileup", help="count buffy support at somatic sites")
    r.add_argument("--vcf", required=True)
    r.add_argument("--bam", required=True)
    r.add_argument("--ref", required=False)       # accepted; pysam pileup needs none
    r.add_argument("--min-bq", dest="min_bq", type=int, default=20)
    r.add_argument("--min-mapq", dest="min_mapq", type=int, default=20)
    r.add_argument("--out", required=True)
    r.set_defaults(func=cmd_pileup)

    s = sub.add_parser("selftest")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
