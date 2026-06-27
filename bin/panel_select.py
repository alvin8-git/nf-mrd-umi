#!/usr/bin/env python3
"""panel_select.py — Pipeline A / Stage A6: personalized panel selection.

Turn a patient's annotated somatic WES variants into the MRD panel (BED + VCF)
that Pipeline B interrogates. The selection rules ARE the panel-design IP:

  1. SNVs only (indel interrogation is harder/noisier — flagged, not selected).
  2. Drop CHIP-associated variants (design P4): gene on the CHIP blocklist OR
     present in the matched buffy-coat normal. CHIP is the #1 cfDNA false
     positive and germline filtering alone does not catch it.
  3. Prefer truncal, high-CCF, clonal variants (design P6): rank by clonal
     probability then CCF, so the panel tracks the tumor trunk, not branches.
  4. Cap at panel_size, enforce a min CCF floor.

Run       : panel_select.py run --vcf somatic.vep.vcf --ccf pyclone.tsv \
              --normal-evidence buffy.tsv --chip-blocklist chip_genes.txt \
              --panel-size 50 --out-bed panel.bed --out-vcf panel.vcf
Self-test : panel_select.py selftest
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    ccf: float                 # cancer cell fraction (PyClone-vi)
    clonal_prob: float         # P(clonal / truncal) (PyClone-vi)
    normal_alt: int            # alt molecules in matched buffy-coat normal
    normal_depth: int          # depth in buffy-coat normal
    pop_af: float = 0.0         # gnomAD/dbSNP population allele freq (germline excl)
    probe_designable: bool = True   # allele-specific probe buildable (probe wf)
    probe_specificity: float = 1.0  # predicted probe on-target specificity [0,1]
    predicted_enrichment: float = 1.0  # probe-design predicted enrichment factor

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1

    @property
    def normal_vaf(self) -> float:
        return self.normal_alt / self.normal_depth if self.normal_depth else 0.0


# ---------------------------------------------------------------------------
# Pure selection logic (unit-testable)
# ---------------------------------------------------------------------------
def is_chip(v: Variant, blocklist: set[str]) -> bool:
    return v.gene.upper() in blocklist


def present_in_normal(v: Variant, min_alt: int = 3, min_vaf: float = 0.02) -> bool:
    """Variant has real support in the buffy-coat normal => germline or CHIP.
    Both an absolute molecule count AND a VAF floor must clear, so a couple of
    stray error reads at high depth don't wrongly reject a true somatic site."""
    return v.normal_alt >= min_alt and v.normal_vaf >= min_vaf


def select_variants(variants: list[Variant], blocklist: set[str],
                    panel_size: int = 50, min_ccf: float = 0.10,
                    normal_min_alt: int = 3, normal_min_vaf: float = 0.02,
                    pop_af_max: float = 1e-3, min_probe_specificity: float = 0.90
                    ) -> tuple[list[Variant], dict[str, int]]:
    """Return (selected, reject_counts). Deterministic ordering.

    A variant must be (a) an SNV, (b) not CHIP, (c) absent from buffy, (d) not a
    population germline variant (gnomAD/dbSNP), (e) clonal enough, AND (f)
    carry a designable, specific allele-specific probe (2Strands enrichment).
    The final panel is the intersection of 'informative' and 'enrichable'.
    """
    reasons = {"non_snv": 0, "chip_gene": 0, "in_normal": 0, "germline_pop": 0,
               "low_ccf": 0, "not_probeable": 0, "probe_offtarget": 0}
    candidates: list[Variant] = []
    for v in variants:
        if not v.is_snv:
            reasons["non_snv"] += 1
            continue
        if is_chip(v, blocklist):
            reasons["chip_gene"] += 1
            continue
        if present_in_normal(v, normal_min_alt, normal_min_vaf):
            reasons["in_normal"] += 1
            continue
        if v.pop_af >= pop_af_max:                  # gnomAD/dbSNP germline exclusion
            reasons["germline_pop"] += 1
            continue
        if v.ccf < min_ccf:
            reasons["low_ccf"] += 1
            continue
        if not v.probe_designable:                  # external probe-design verdict
            reasons["not_probeable"] += 1
            continue
        if v.probe_specificity < min_probe_specificity:
            reasons["probe_offtarget"] += 1
            continue
        candidates.append(v)
    # Truncal AND most-enrichable first: clonal_prob, CCF, predicted enrichment.
    candidates.sort(key=lambda v: (-v.clonal_prob, -v.ccf,
                                   -v.predicted_enrichment, v.chrom, v.pos))
    return candidates[:panel_size], reasons


def probe_prescreen(context_seq: str, var_offset: int, ref: str, alt: str,
                    probe_len: int = 30) -> dict:
    """Heuristic allele-specific probe feasibility, run BEFORE the external
    probe-design workflow to pre-rank candidates. ponytail: GC + homopolymer + N
    proxies only — real Tm, secondary structure, and off-target are the existing
    probe-design tool's job; this just avoids shipping obviously-undesignable
    sites into it. Returns {designable, gc, max_homopolymer, reason}.
    """
    seq = context_seq.upper()
    half = probe_len // 2
    window = seq[max(0, var_offset - half): var_offset + half + 1]
    if not window or "N" in window:
        return {"designable": False, "gc": 0.0, "max_homopolymer": 0,
                "reason": "N_or_empty"}
    gc = (window.count("G") + window.count("C")) / len(window)
    max_run = run = 1
    for i in range(1, len(window)):
        run = run + 1 if window[i] == window[i - 1] else 1
        max_run = max(max_run, run)
    gc_ok = 0.30 <= gc <= 0.70
    designable = gc_ok and max_run < 5
    reason = "ok" if designable else ("gc" if not gc_ok else "homopolymer")
    return {"designable": designable, "gc": round(gc, 3),
            "max_homopolymer": max_run, "reason": reason}


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def load_variants(vcf_path: str, ccf_tsv: str, normal_tsv: str) -> list[Variant]:
    """Join somatic VCF (gene from VEP CSQ) with PyClone-vi CCF and buffy
    normal evidence, all keyed on (chrom,pos,ref,alt)."""
    import pysam

    ccf = {_k(r): (float(r["ccf"]), float(r.get("clonal_prob", 0.0)))
           for r in _read_tsv(ccf_tsv)}
    norm = {_k(r): (int(r["normal_alt"]), int(r["normal_depth"]))
            for r in _read_tsv(normal_tsv)}

    out: list[Variant] = []
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            for alt in rec.alts or []:
                key = (rec.chrom, rec.pos, rec.ref, alt)
                gene = _gene_from_csq(rec)
                c, cp = ccf.get(key, (0.0, 0.0))
                na, nd = norm.get(key, (0, 0))
                out.append(Variant(rec.chrom, rec.pos, rec.ref, alt, gene,
                                   c, cp, na, nd))
    return out


def _k(r: dict) -> tuple[str, int, str, str]:
    return (r["chrom"], int(r["pos"]), r["ref"], r["alt"])


def _read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _gene_from_csq(rec) -> str:
    """Pull SYMBOL from the first VEP CSQ entry. ponytail: positional parse of
    the first transcript is enough for gene-level CHIP blocklisting."""
    csq = rec.info.get("CSQ")
    if not csq:
        return ""
    first = (csq[0] if isinstance(csq, (list, tuple)) else csq).split("|")
    # VEP default field order: Allele|Consequence|IMPACT|SYMBOL|...
    return first[3] if len(first) > 3 else ""


def write_bed(selected: list[Variant], path: str) -> None:
    with open(path, "w") as fh:
        for v in selected:
            name = f"{v.gene}:{v.ref}>{v.alt}"
            fh.write(f"{v.chrom}\t{v.pos - 1}\t{v.pos}\t{name}\n")  # 0-based BED


def write_vcf(selected: list[Variant], src_vcf: str, path: str) -> None:
    """Subset the source VCF to the selected sites, preserving its header."""
    import pysam

    keep = {(v.chrom, v.pos, v.ref, v.alt) for v in selected}
    with pysam.VariantFile(src_vcf) as src:
        with pysam.VariantFile(path, "w", header=src.header) as out:
            for rec in src:
                for alt in rec.alts or []:
                    if (rec.chrom, rec.pos, rec.ref, alt) in keep:
                        out.write(rec)
                        break


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def cmd_run(args) -> None:
    blocklist = {g.strip().upper() for g in open(args.chip_blocklist)
                 if g.strip() and not g.startswith("#")}
    variants = load_variants(args.vcf, args.ccf, args.normal_evidence)
    selected, reasons = select_variants(
        variants, blocklist, panel_size=args.panel_size, min_ccf=args.min_ccf)
    if not selected:
        sys.exit("ERROR: 0 variants selected — panel would be empty. "
                 f"Inspect rejects: {reasons}")
    write_bed(selected, args.out_bed)
    write_vcf(selected, args.vcf, args.out_vcf)
    print(f"panel_select: {len(selected)}/{len(variants)} selected "
          f"-> {args.out_bed}, {args.out_vcf}  rejects={reasons}")


def cmd_selftest(_args) -> None:
    blocklist = {"DNMT3A", "TET2", "ASXL1"}
    V = Variant
    variants = [
        V("chr1", 100, "C", "T", "TP53",   0.95, 0.98, 0, 4000),  # truncal -> keep
        V("chr1", 200, "G", "A", "DNMT3A", 0.90, 0.95, 0, 4000),  # CHIP gene -> drop
        V("chr2", 300, "A", "G", "KRAS",   0.40, 0.30, 0, 4000),  # subclonal -> keep (low rank)
        V("chr3", 400, "C", "T", "EGFR",   0.92, 0.97, 120, 4000),# in normal -> drop
        V("chr4", 500, "AT", "A", "APC",   0.99, 0.99, 0, 4000),  # indel -> drop
        V("chr5", 600, "C", "T", "PIK3CA", 0.05, 0.50, 0, 4000),  # low CCF -> drop
        V("chr6", 700, "C", "T", "BRCA1",  0.90, 0.90, 0, 4000,
          pop_af=0.05),                                           # germline pop -> drop
        V("chr7", 800, "C", "T", "BRAF",   0.90, 0.90, 0, 4000,
          probe_designable=False),                               # un-probeable -> drop
    ]
    selected, reasons = select_variants(variants, blocklist, panel_size=50,
                                        min_ccf=0.10)
    genes = [v.gene for v in selected]
    assert genes == ["TP53", "KRAS"], genes          # everything non-informative/un-probeable gone
    assert selected[0].gene == "TP53", "truncal must rank first"
    assert reasons == {"non_snv": 1, "chip_gene": 1, "in_normal": 1,
                       "germline_pop": 1, "low_ccf": 1,
                       "not_probeable": 1, "probe_offtarget": 0}, reasons

    # panel_size cap respected.
    capped, _ = select_variants(variants, blocklist, panel_size=1, min_ccf=0.10)
    assert capped == [variants[0]], capped

    # present_in_normal needs BOTH count and vaf.
    v_lowvaf = V("chrX", 1, "C", "T", "X", 0.9, 0.9, 5, 100000)   # 5 alt, vaf 5e-5
    assert not present_in_normal(v_lowvaf), "high-depth noise should not reject"

    # Predicted enrichment breaks ties: more-enrichable probe ranks first.
    e_lo = V("chr8", 1, "C", "T", "LO", 0.9, 0.9, 0, 4000, predicted_enrichment=10.0)
    e_hi = V("chr8", 2, "C", "T", "HI", 0.9, 0.9, 0, 4000, predicted_enrichment=80.0)
    sel_e, _ = select_variants([e_lo, e_hi], blocklist, panel_size=2, min_ccf=0.10)
    assert [v.gene for v in sel_e] == ["HI", "LO"], sel_e

    # Probe pre-screen: balanced sequence designable; homopolymer / N not.
    assert probe_prescreen("ACGT" * 10, 20, "C", "T")["designable"] is True
    assert probe_prescreen("A" * 40, 20, "C", "T")["designable"] is False
    assert probe_prescreen("ACGTNACGT", 4, "C", "T", probe_len=6)["designable"] is False

    print("panel_select.py selftest: PASS "
          f"(kept {genes}, rejects {reasons})")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="select the personalized panel")
    r.add_argument("--vcf", required=True, help="annotated somatic VCF (VEP)")
    r.add_argument("--ccf", required=True, help="PyClone-vi CCF TSV")
    r.add_argument("--normal-evidence", dest="normal_evidence", required=True,
                   help="buffy-coat alt/depth TSV")
    r.add_argument("--chip-blocklist", dest="chip_blocklist", required=True)
    r.add_argument("--panel-size", dest="panel_size", type=int, default=50)
    r.add_argument("--min-ccf", dest="min_ccf", type=float, default=0.10)
    r.add_argument("--out-bed", dest="out_bed", required=True)
    r.add_argument("--out-vcf", dest="out_vcf", required=True)
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("selftest", help="run selection-logic asserts")
    s.set_defaults(func=cmd_selftest)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
