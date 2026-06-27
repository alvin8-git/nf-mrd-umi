# Validation Runbook

How to take this pipeline out of the silo: run it against real public data
(SEQC2 liquid biopsy, BioProject PRJNA677999) and get real specificity, LoD, and
LoQ numbers instead of simulated ones. Every dataset and accession here was
verified against NCBI/EGA (see [VALIDATION_DATASETS.md](VALIDATION_DATASETS.md)).

## What this validates (and what it does not)

Validates: UMI/consensus, interrogation, the empirical null, and the LoB/LoD/LoQ
statistics, on real contrived ctDNA at known VAFs.

Does NOT validate: the allele-specific enrichment layer. SEQC2 uses hybrid/
amplicon assays, not 2Strands probes, so score with enrichment = 1. The
enrichment de-biasing and per-probe calibration still need in-house dilutions
through the real probes. This runbook covers everything else.

## The chain

```
PRJNA677999 (SEQC2)                              verified on NCBI
  Step 1  build_panel.sh        -> panel.vcf.gz  truth SNVs INTERSECT assay BED
  Step 2  fetch_and_interrogate -> results/interrogate/<RUN>.tsv  per-run counts
  Step 3  build background+PoN  -> background.tsv, pon.tsv  from blank runs
  Step 4  validate.py run-real  -> results/seqc2_validation.json  FPR/LoD/LoQ
```

Cohort (TFS2 assay, 3 labs), from
[data/manifests/seqc2_tfs2_titration.manifest.tsv](data/manifests/seqc2_tfs2_titration.manifest.tsv):
12 blanks (VAF 0), 11 at 0.4%, 24 at 2%, 12 at 10%.

---

## Prerequisites

```
sra-tools (prefetch, fasterq-dump), bwa-mem2, samtools, bcftools, bedtools, tabix
python3 with numpy, scipy, pysam
fgbio (+ Java)   # only if running the UMI duplex path
```

You also supply:
- **GRCh38 reference**, bwa-mem2-indexed (`bwa-mem2 index GRCh38.fa`). Match the
  build the SEQC2 truth set uses (GRCh38); harmonize chrom naming in Step 1 if needed.
- **TFS2 assay target BED** - the regions the assay covers. Source: the SEQC2
  Liquid Biopsy figshare project (per-vendor panel region BEDs). This is the one
  input with no auto-download.

Disk: budget ~tens of GB for the 59 runs' intermediates (cleaned automatically
unless you pass `-k`).

---

## Step 1 - Build the tumor-informed panel

The "panel" is the known HCC1395 mutations the assay can see: SEQC2 high-confidence
somatic SNVs, intersected with the confident regions and the assay BED.

```bash
bin/build_panel.sh -b TFS2_targets.bed -o panel -r GRCh38.fa
# add -C ucsc (or -C ensembl) if chrom naming differs from your reference
```

Outputs `panel.vcf.gz` (+ `panel.bed`). It auto-downloads the verified SEQC2
truth set. **Check the site count it prints** - if it warns that the panel is
small (<8 sites), this hotspot assay may not cover enough HCC1395 mutations for
panel-integrated MRD; consider a broader assay's BED.

## Step 2 - Download and interrogate every run

```bash
bin/fetch_and_interrogate.sh \
  -m data/manifests/seqc2_tfs2_titration.manifest.tsv \
  -r GRCh38.fa -p panel.vcf.gz -j 4 -t 8
# UMI duplex path (if the assay carries UMIs): add -u -s "8M+T 8M+T"
#   (read structure is assay-specific; do not guess it)
```

Fills `results/interrogate/<RUN>.tsv` for all 59 runs. Resumable. Without `-u`,
duplicates are flagged but not collapsed, so molecule counts are inflated; for a
UMI assay use `-u` for true molecule-level counts.

## Step 3 - Build the background and PoN from blank runs

The empirical null and the per-site error model come from tumor-free (Bf, VAF 0)
runs. Keep the folds DISJOINT: fit on half the blanks, test specificity on the
other half.

```bash
# split the 12 blanks: 6 to fit the background, 6 stay as held-out test
grep -P '\tblank\t' data/manifests/seqc2_tfs2_titration.manifest.tsv \
  | cut -f1 > all_blanks.txt
head -6 all_blanks.txt > train_blanks.txt   # background-fit set
tail -6 all_blanks.txt > test_blanks.txt     # held-out specificity set

# drop the training blanks from the run-real manifest (so test stays disjoint)
grep -vF -f train_blanks.txt \
  data/manifests/seqc2_tfs2_titration.manifest.tsv > manifest.test.tsv
```

Then build the background and the empirical-null PoN from the TRAIN blanks'
interrogate outputs:

```bash
python3 bin/build_background.py run \
  --out-background background.tsv --out-pon pon.tsv \
  $(sed 's#^#results/interrogate/#; s#$#.tsv#' train_blanks.txt)
```

Outputs:
- **`background.tsv`** (`chrom pos ref alt alpha beta enrich`) - per-site
  beta-binomial fit by method-of-moments over the train blanks; `enrich = 1`
  (TFS2 is not enrichment-based).
- **`pon.tsv`** (`chrom pos ref alt <blank1> <blank2> ...`) - per-donor error
  rates, the covariance-preserving empirical null.

## Step 4 - Score the cohort

```bash
python3 bin/validate.py run-real \
  --manifest manifest.test.tsv \
  --background background.tsv \
  --pon pon.tsv \
  --out results/seqc2_validation.json
```

## Step 5 - Read the result

The report gives the first REAL numbers:

- **observed FPR** on the held-out blanks vs nominal alpha (0.05). Near alpha =
  the empirical null holds on real data. Much higher = the null is still
  anticonservative on this assay's real error structure.
- **LoD95** - lowest VAF detected in >=95% of replicates. Compare to the 0.4%
  (Ff) level, the lowest in this cohort.
- **LoQ** + per-level `rel_bias` - how well tumor fraction tracks the known VAF.
- **hit_rate** per level - sensitivity at 0.4% / 2% / 10%.

Compare against `bin/validate.py simulate-run` (the simulated baseline). Real
numbers worse than simulated is the expected, honest outcome: real data has
artifacts the simulator does not model. That gap is the next round of work.

---

## Caveats (carried from the data, not bugs)

1. **Enrichment layer is not tested here** (no allele-specific probes in public
   data). Score with enrichment = 1.
2. **`vaf_truth` is approximate** (`0.5 x mixture fraction`, het assumption). Fine
   for LoD binning; recompute from the truth-set VAF distribution over your panel
   for precise quantification bias.
3. **Input mass varies** (Ff is 100 ng, others 10-25 ng) - a depth covariate;
   subset to one input for the cleanest read.
4. **No-UMI molecule counting** inflates depth (duplicates flagged, not
   collapsed). Use `-u` for UMI assays.
5. **Cohort hygiene** - never fit the background and test specificity on the same
   blanks (Step 3 enforces this).

## See also

- [VALIDATION_DATASETS.md](VALIDATION_DATASETS.md) - full dataset catalogue
- [data/manifests/README.md](data/manifests/README.md) - the cohort and its caveats
- [TODO.md](TODO.md) - P0 "validation harness on real data" and the open items
- [Documentation.md](Documentation.md) - pipeline 101
