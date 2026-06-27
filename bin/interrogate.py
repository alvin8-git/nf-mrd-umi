#!/usr/bin/env python3
"""interrogate.py — Pipeline B / Stage B5: known-site interrogation.

NOT variant calling (design P1). For each patient-specific panel site, count how
many UNIQUE MOLECULES (consensus reads — one consensus read == one source
molecule) support the tumor's known ALT allele above background. Emits the
per-site table that mrd_integrate.py consumes.

Input  : molecularly-collapsed, realigned consensus BAM (Stage B4) + panel VCF.
Output : TSV with one row per site:
         chrom pos ref alt depth_unique_molecules alt_molecule_count
         duplex_support mean_consensus_qual

Run       : interrogate.py run --bam cons.bam --panel panel.vcf --out counts.tsv
Self-test : interrogate.py selftest   (tests the pure pileup-counting core;
            no BAM required)
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, asdict

_ACGT = frozenset("ACGT")


@dataclass
class SiteCount:
    chrom: str
    pos: int
    ref: str
    alt: str
    depth_unique_molecules: int
    alt_molecule_count: int
    duplex_support: int        # alt molecules backed by BOTH strands (duplex)
    mean_consensus_qual: float


# ---------------------------------------------------------------------------
# Pure core (unit-testable without pysam / a BAM)
# ---------------------------------------------------------------------------
def count_site(bases, quals, is_duplex, ref: str, alt: str,
               min_bq: int = 30) -> SiteCount:
    """Tally one pileup column.

    bases     : per-molecule called base at the site (str, uppercased)
    quals     : per-molecule consensus base quality (int)
    is_duplex : per-molecule bool — True if both source strands agreed (fgbio
                duplex consensus: per-strand depths aD>0 and bD>0)
    Counts only A/C/G/T bases passing min_bq. depth = qualifying molecules;
    alt_count = those equal to ALT; duplex_support = ALT molecules that are
    duplex. mean_consensus_qual averaged over qualifying molecules.
    """
    if not (len(bases) == len(quals) == len(is_duplex)):
        raise ValueError("bases/quals/is_duplex length mismatch")
    alt = alt.upper()
    depth = alt_count = duplex_support = 0
    qual_sum = 0
    for base, q, dup in zip(bases, quals, is_duplex):
        b = base.upper()
        if b not in _ACGT or q < min_bq:
            continue
        depth += 1
        qual_sum += q
        if b == alt:
            alt_count += 1
            if dup:
                duplex_support += 1
    mean_q = (qual_sum / depth) if depth else 0.0
    return SiteCount("", 0, ref, alt, depth, alt_count, duplex_support,
                     round(mean_q, 1))


# ---------------------------------------------------------------------------
# pysam wrapper (thin; not exercised by selftest)
# ---------------------------------------------------------------------------
def _read_duplex_flag(read) -> bool:
    """fgbio CallDuplexConsensusReads sets per-strand raw depths aD/bD.
    A molecule is duplex when both strands contributed (aD>0 and bD>0)."""
    try:
        return read.get_tag("aD") > 0 and read.get_tag("bD") > 0
    except KeyError:
        return False


def interrogate_bam(bam_path: str, panel_vcf: str, out_tsv: str,
                    min_bq: int = 30, min_mapq: int = 20) -> None:
    import pysam  # ponytail: imported here so selftest needs no pysam install

    sites: list[tuple[str, int, str, str]] = []
    with pysam.VariantFile(panel_vcf) as vcf:
        for rec in vcf:
            for a in rec.alts or []:
                sites.append((rec.chrom, rec.pos, rec.ref, a))

    bam = pysam.AlignmentFile(bam_path, "rb")
    rows: list[SiteCount] = []
    for chrom, pos1, ref, alt in sites:
        pos0 = pos1 - 1  # VCF 1-based -> pysam 0-based
        bases, quals, dup = [], [], []
        for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True,
                              min_base_quality=0, stepper="samtools"):
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.mapping_quality < min_mapq:
                    continue
                bases.append(aln.query_sequence[pr.query_position])
                quals.append(aln.query_qualities[pr.query_position])
                dup.append(_read_duplex_flag(aln))
        sc = count_site(bases, quals, dup, ref, alt, min_bq=min_bq)
        sc.chrom, sc.pos = chrom, pos1
        rows.append(sc)
    bam.close()

    _write_tsv(rows, out_tsv)
    print(f"interrogate: {len(rows)} sites -> {out_tsv}")


def _write_tsv(rows: list[SiteCount], out_tsv: str) -> None:
    fields = ["chrom", "pos", "ref", "alt", "depth_unique_molecules",
              "alt_molecule_count", "duplex_support", "mean_consensus_qual"]
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def cmd_run(args) -> None:
    interrogate_bam(args.bam, args.panel, args.out,
                    min_bq=args.min_bq, min_mapq=args.min_mapq)


def cmd_selftest(_args) -> None:
    # 10 molecules: 6 ref(C), 3 alt(T), 1 low-qual alt(dropped). 2 alt are duplex.
    bases = ["C", "C", "C", "C", "C", "C", "T", "T", "T", "T"]
    quals = [35,  35,  35,  35,  35,  35,  35,  35,  35,  10]   # last dropped
    dup =   [False]*6 + [True, True, False, False]
    sc = count_site(bases, quals, dup, ref="C", alt="T", min_bq=30)
    assert sc.depth_unique_molecules == 9, sc          # 10 - 1 low-qual
    assert sc.alt_molecule_count == 3, sc              # 3 high-qual T
    assert sc.duplex_support == 2, sc
    assert abs(sc.mean_consensus_qual - 35.0) < 1e-6, sc

    # All ref -> zero alt.
    sc0 = count_site(["C"]*8, [40]*8, [True]*8, "C", "T")
    assert sc0.alt_molecule_count == 0 and sc0.depth_unique_molecules == 8, sc0

    # Empty column -> zero depth, zero mean qual, no crash.
    sce = count_site([], [], [], "C", "T")
    assert sce.depth_unique_molecules == 0 and sce.mean_consensus_qual == 0.0, sce

    # Non-ACGT (N) ignored.
    scn = count_site(["N", "T"], [40, 40], [False, False], "C", "T")
    assert scn.depth_unique_molecules == 1 and scn.alt_molecule_count == 1, scn

    # Length-mismatch guard.
    try:
        count_site(["C"], [40, 40], [True], "C", "T")
        raise AssertionError("expected ValueError on length mismatch")
    except ValueError:
        pass

    print("interrogate.py selftest: PASS")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="count molecules at panel sites")
    r.add_argument("--bam", required=True, help="consensus BAM (Stage B4)")
    r.add_argument("--panel", required=True, help="panel VCF (ref/alt per site)")
    r.add_argument("--out", required=True)
    r.add_argument("--min-bq", dest="min_bq", type=int, default=30)
    r.add_argument("--min-mapq", dest="min_mapq", type=int, default=20)
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("selftest", help="test the pure counting core")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
