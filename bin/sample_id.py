#!/usr/bin/env python3
"""sample_id.py - sample-identity + patient-lock for tumor-informed MRD.

Closes the P0 CAP/CLIA gap (review-findings): prove the cfDNA, tumor, and buffy
BAMs are the SAME patient, and that a cfDNA run is interrogated against ITS OWN
patient's panel (not a swapped one).

Two mechanisms:
  1. SNP fingerprint concordance - genotype a fixed common-SNP set in each BAM,
     compare genotype dosage across samples. Low concordance => sample swap.
  2. Patient-lock provenance hash - a deterministic token over the panel's sites
     (+ patient id). Stamp it at panel build; verify it before interrogation so a
     cfDNA sample can never be scored against another patient's panel.

Subcommands:
  fingerprint --bam --sites --out [--min-depth] [--min-mapq] [--min-bq]
  concordance --a --b [--out] [--same-thresh] [--diff-thresh]
  provenance  --panel [--patient-id] [--out]      # emit lock token
  verify-lock --panel --lock                       # assert panel matches token
  selftest

The fingerprint SNP set is assay-independent: provide a standard common-SNP list
(e.g. a ~100-SNP fingerprinting panel) as a VCF/BED via --sites. We do NOT bundle
one; document the source in VALIDATION.md.
"""
import argparse
import csv
import hashlib
import sys


# ---------- pure cores (selftest-covered; no I/O) ----------

def dosage(alt_frac, min_het=0.15, max_het=0.85):
    """Alt-allele fraction -> diploid genotype dosage {0,1,2} or None (no-call).

    ponytail: fixed-threshold genotyper, not a likelihood model. Fingerprinting
    only needs to tell hom-ref / het / hom-alt apart at common SNPs; a GL model
    is upgrade-path if borderline calls ever drive a false swap flag.
    """
    if alt_frac is None:
        return None
    if alt_frac < min_het:
        return 0
    if alt_frac > max_het:
        return 2
    return 1


def concordance(geno_a, geno_b):
    """Two {site: dosage} maps -> (concordant_fraction, n_compared).

    Compares only sites genotyped (non-None) in BOTH. Identical-twin / same-sample
    pairs ~1.0; unrelated ~0.4-0.6 (chance agreement of three genotype classes).
    """
    common = [s for s in geno_a if s in geno_b
              and geno_a[s] is not None and geno_b[s] is not None]
    if not common:
        return (None, 0)
    match = sum(1 for s in common if geno_a[s] == geno_b[s])
    return (match / len(common), len(common))


def verdict(frac, n, same_thresh=0.90, diff_thresh=0.70, min_sites=20):
    """Concordance -> SAME / DIFFERENT / AMBIGUOUS (too few sites or middling)."""
    if frac is None or n < min_sites:
        return "AMBIGUOUS"
    if frac >= same_thresh:
        return "SAME"
    if frac < diff_thresh:
        return "DIFFERENT"
    return "AMBIGUOUS"


def panel_token(sites, patient_id=""):
    """Deterministic patient-lock hash over sorted (chrom,pos,ref,alt) + patient.

    Order-independent (sites are sorted first) so the token depends only on the
    panel CONTENT, not VCF line order. Same panel + same patient => same token.
    """
    norm = sorted("{}:{}:{}:{}".format(c, p, r, a) for c, p, r, a in sites)
    h = hashlib.sha256()
    h.update(patient_id.encode())
    h.update(b"\x00")
    h.update("\n".join(norm).encode())
    return "panel-{}-{}".format(len(norm), h.hexdigest()[:16])


# ---------- I/O wrappers (need pysam / files; not in selftest) ----------

def read_sites(path):
    """VCF(.gz) or BED-ish TSV -> list of (chrom, pos1, ref, alt). pos is 1-based."""
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            # VCF: chrom pos id ref alt ...   BED4: chrom pos ref alt (pos 1-based)
            if len(f) >= 5 and f[1].isdigit():
                out.append((f[0], int(f[1]), f[3], f[4]))
            elif len(f) >= 4 and f[1].isdigit():
                out.append((f[0], int(f[1]), f[2], f[3]))
    return out


def fingerprint_bam(bam, sites, min_depth=10, min_mapq=20, min_bq=20):
    """Genotype each SNP site in a BAM -> {site_key: dosage}. Needs pysam."""
    import pysam
    sf = pysam.AlignmentFile(bam, "rb")
    geno, rows = {}, []
    for chrom, pos, ref, alt in sites:
        key = "{}:{}".format(chrom, pos)
        depth = alt_n = 0
        for col in sf.pileup(chrom, pos - 1, pos, truncate=True,
                             min_base_quality=min_bq, min_mapping_quality=min_mapq):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                b = pr.alignment.query_sequence[pr.query_position].upper()
                depth += 1
                if b == alt.upper():
                    alt_n += 1
        frac = (alt_n / depth) if depth >= min_depth else None
        g = dosage(frac)
        geno[key] = g
        rows.append((chrom, pos, ref, alt, depth, alt_n,
                     "" if frac is None else "{:.3f}".format(frac),
                     "" if g is None else g))
    sf.close()
    return geno, rows


def load_fingerprint(path):
    """Read a fingerprint TSV back into {site_key: dosage}."""
    geno = {}
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r, None)
        for row in r:
            if not row:
                continue
            key = "{}:{}".format(row[0], row[1])
            g = row[7]
            geno[key] = int(g) if g not in ("", "None") else None
    return geno


# ---------- subcommands ----------

def cmd_fingerprint(a):
    sites = read_sites(a.sites)
    if not sites:
        sys.exit("ERROR: no sites read from {}".format(a.sites))
    _, rows = fingerprint_bam(a.bam, sites, a.min_depth, a.min_mapq, a.min_bq)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "depth", "alt", "alt_frac", "geno"])
        w.writerows(rows)
    called = sum(1 for r in rows if r[7] != "")
    print("[fingerprint] {} sites, {} genotyped -> {}".format(len(rows), called, a.out))


def cmd_concordance(a):
    ga, gb = load_fingerprint(a.a), load_fingerprint(a.b)
    frac, n = concordance(ga, gb)
    v = verdict(frac, n, a.same_thresh, a.diff_thresh)
    line = "concordance={} sites={} verdict={}".format(
        "NA" if frac is None else "{:.4f}".format(frac), n, v)
    print(line)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(line + "\n")
    # non-zero exit on a confirmed swap so a pipeline step can gate on it
    sys.exit(1 if v == "DIFFERENT" else 0)


def cmd_provenance(a):
    sites = read_sites(a.panel)
    tok = panel_token(sites, a.patient_id)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(tok + "\n")
    print(tok)


def cmd_verify_lock(a):
    sites = read_sites(a.panel)
    with open(a.lock) as fh:
        want = fh.read().strip()
    # Token is opaque (no patient id recoverable from it), so the caller passes the
    # same --patient-id used at lock time; recompute over the panel and compare.
    got = panel_token(sites, a.patient_id)
    if got != want:
        sys.exit("LOCK MISMATCH: panel does not match patient lock\n  want {}\n  got  {}"
                 .format(want, got))
    print("LOCK OK: {}".format(got))


def selftest(_a):
    # dosage thresholds
    assert dosage(0.0) == 0 and dosage(0.5) == 1 and dosage(1.0) == 2
    assert dosage(0.14) == 0 and dosage(0.86) == 2 and dosage(None) is None
    # concordance: identical -> 1.0; disjoint coverage -> only co-covered count
    same = {"1:100": 0, "1:200": 1, "1:300": 2}
    assert concordance(same, dict(same)) == (1.0, 3)
    half = {"1:100": 0, "1:200": 1, "1:300": 0}  # one mismatch of three
    f, n = concordance(same, half)
    assert n == 3 and abs(f - 2/3) < 1e-9
    # no-calls excluded from the denominator
    withN = {"1:100": 0, "1:200": None, "1:300": 2}
    assert concordance(same, withN) == (1.0, 2)
    # verdicts
    assert verdict(0.99, 50) == "SAME"
    assert verdict(0.50, 50) == "DIFFERENT"
    assert verdict(0.80, 50) == "AMBIGUOUS"
    assert verdict(1.0, 5) == "AMBIGUOUS"      # too few sites
    assert verdict(None, 0) == "AMBIGUOUS"
    # patient-lock: order-independent, patient- and content-sensitive
    s1 = [("1", 100, "A", "G"), ("2", 200, "C", "T")]
    s2 = [("2", 200, "C", "T"), ("1", 100, "A", "G")]  # reordered
    assert panel_token(s1, "P1") == panel_token(s2, "P1")     # order-independent
    assert panel_token(s1, "P1") != panel_token(s1, "P2")     # patient-sensitive
    assert panel_token(s1, "P1") != panel_token(s1[:1], "P1")  # content-sensitive
    print("sample_id selftest OK")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fingerprint", help="genotype a common-SNP set in a BAM")
    f.add_argument("--bam", required=True)
    f.add_argument("--sites", required=True, help="common-SNP VCF(.gz) or BED4")
    f.add_argument("--out", required=True)
    f.add_argument("--min-depth", dest="min_depth", type=int, default=10)
    f.add_argument("--min-mapq", dest="min_mapq", type=int, default=20)
    f.add_argument("--min-bq", dest="min_bq", type=int, default=20)
    f.set_defaults(func=cmd_fingerprint)

    c = sub.add_parser("concordance", help="compare two fingerprints; exit 1 if DIFFERENT")
    c.add_argument("--a", required=True)
    c.add_argument("--b", required=True)
    c.add_argument("--out")
    c.add_argument("--same-thresh", dest="same_thresh", type=float, default=0.90)
    c.add_argument("--diff-thresh", dest="diff_thresh", type=float, default=0.70)
    c.set_defaults(func=cmd_concordance)

    pr = sub.add_parser("provenance", help="emit patient-lock token for a panel")
    pr.add_argument("--panel", required=True)
    pr.add_argument("--patient-id", dest="patient_id", default="")
    pr.add_argument("--out")
    pr.set_defaults(func=cmd_provenance)

    v = sub.add_parser("verify-lock", help="assert panel matches a lock token (exit!=0 on mismatch)")
    v.add_argument("--panel", required=True)
    v.add_argument("--lock", required=True)
    v.add_argument("--patient-id", dest="patient_id", default="")
    v.set_defaults(func=cmd_verify_lock)

    s = sub.add_parser("selftest", help="assert pure-core mechanics")
    s.set_defaults(func=selftest)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
