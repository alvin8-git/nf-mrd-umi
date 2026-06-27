# run-real manifest: SEQC2 TFS2 titration (PRJNA677999)

`seqc2_tfs2_titration.manifest.tsv` is a real-data validation cohort for
`validate.py run-real`, drawn from the SEQC2 Liquid Biopsy study
(BioProject PRJNA677999, NCBI-verified). It is a tumor-into-normal titration of
the HCC1395 (tumor "Sample A") / HCC1395BL (normal "Sample B") cell-line pair,
sequenced by one ctDNA assay (TFS2) across three labs.

## What's in it (59 runs)

| sample_code | mixture | vaf_truth | label | runs |
|---|---|---|---|---|
| Bf | 100% normal | 0 | blank | 12 |
| Ff | 0.8% tumor | 0.004 | tumor | 11 |
| Ef | 4% tumor | 0.02 | tumor | 24 |
| Df | 20% tumor | 0.10 | tumor | 12 |

Columns: `sample` (SRA run), `label` (blank/tumor), `vaf_truth`,
`site_counts_path` (where `interrogate.py` output will land), plus provenance
(`sample_code`, `mixture_fraction_A`, `assay`, `lab`, `input`, `biosample`).
`validate.py run-real` reads the first four; the rest document the cohort.

## Important caveats (read before trusting any number)

1. **`vaf_truth` is the expected panel-mean VAF, approximated as
   `0.5 x mixture_fraction_A`** (heterozygous-variant assumption, since most
   HCC1395 somatic SNVs are heterozygous). The exact per-variant truth is the
   SEQC2 high-confidence HCC1395 call set VAF scaled by `mixture_fraction_A`.
   For LoD/hit-rate binning the approximation is fine; for precise
   quantification-bias (`rel_bias`), recompute `vaf_truth` per sample as the
   mean truth-set VAF over your panel variants times the fraction.
2. **Input mass is not constant.** Ff (lowest VAF) is 100 ng while the blanks,
   Ef and Df are 10-25 ng. Different input = different unique-molecule depth,
   which affects LoD. For the cleanest read, subset to a single `input` value,
   or normalize. As-is, treat input as a covariate.
3. **TFS2 is a hybrid-capture/amplicon ctDNA assay, NOT allele-specific
   enrichment.** So score with **enrichment factor = 1** (no de-biasing). This
   cohort validates the consensus, interrogation, empirical null, and LoD/LoQ
   math, but it CANNOT validate the 2Strands enrichment layer. That still needs
   in-house dilutions through the real probes.
4. **Cohort hygiene.** Split the 12 Bf blanks: fit the background/PoN on one
   half, evaluate FPR on the held-out half. Never fit and test on the same
   blanks.

## How to fill it in (the upstream pipeline)

`validate.py run-real` needs each row's `site_counts_path` to exist. Produce them:

```bash
# 0. tooling: sra-tools, bwa-mem2, samtools, (fgbio if the assay carries UMIs)
# 1. build the tumor-informed panel for this cohort:
#    SEQC2 high-confidence HCC1395 somatic VCF  INTERSECT  TFS2 panel BED
#    -> panel.vcf   (the variants to interrogate; truth from PRJNA489865)

# 2. per run: download -> align -> (consensus) -> interrogate
while read -r srr; do
  prefetch "$srr" && fasterq-dump --split-files "$srr"
  bwa-mem2 mem -t8 GRCh38.fa "${srr}_1.fastq" "${srr}_2.fastq" \
    | samtools sort -o "${srr}.bam" && samtools index "${srr}.bam"
  # if UMIs present: fgbio GroupReadsByUmi + CallMolecularConsensusReads first
  python3 bin/interrogate.py run \
    --bam "${srr}.bam" --panel panel.vcf --out "results/interrogate/${srr}.tsv"
done < <(tail -n +2 data/manifests/seqc2_tfs2_titration.manifest.tsv | cut -f1)

# 3. build the background (beta-binomial + PoN) from HALF the Bf blanks
#    -> background.tsv  (+ optional pon.tsv for the empirical null)

# 4. score the cohort
python3 bin/validate.py run-real \
  --manifest data/manifests/seqc2_tfs2_titration.manifest.tsv \
  --background background.tsv \
  --out results/seqc2_tfs2_validation.json
```

Expected readout: false-positive rate on the held-out Bf blanks (specificity),
hit-rate per VAF level, LoD95, and quantification bias. These are the first
**real** numbers for the pipeline, replacing the simulated ones in `validate.py
simulate-run`.

## Provenance

Built from NCBI E-utilities runinfo + BioSample attributes for PRJNA677999
(samples SAMN16786373/74/75 and SAMN17007052). Mixture fractions are quoted
verbatim from the BioSample `isolate` attributes. Download each run with
`prefetch <sample>`; no hard-coded URLs (SRA paths are region/version specific).
