# run-real manifest: SEQC2 ILM2 titration (PRJNA677999)

`seqc2_ilm2_titration.manifest.tsv` is the real-data validation cohort for
`validate.py run-real`, drawn from the SEQC2 Liquid Biopsy study
(BioProject PRJNA677999, NCBI-verified). It is a tumor-into-normal titration of
the HCC1395 (tumor "Sample A") / HCC1395BL (normal "Sample B") cell-line pair,
sequenced by the **ILM2** (Illumina, paired-end) ctDNA assay.

> An earlier `seqc2_tfs2_titration.manifest.tsv` (TFS2, Ion Torrent) was
> **abandoned**: those runs are single-end, so the paired-end `fetch_and_interrogate.sh`
> path produced no usable fastq. ILM2 (paired Illumina) is the active cohort.

## What's in it (48 runs)

| sample_code | mixture | vaf_truth | label | runs |
|---|---|---|---|---|
| Bf | 100% normal | 0 | blank | 12 |
| Ef | 4% tumor | 0.02 | tumor | 24 |
| Df | 20% tumor | 0.10 | tumor | 12 |

Columns: `sample` (SRA run), `label` (blank/tumor), `vaf_truth`,
`site_counts_path` (where `interrogate.py` output lands), plus provenance
(`sample_code`, `mixture_fraction_A`, `assay`, `lab`, `input`, `biosample`).
`validate.py run-real` reads the first four; the rest document the cohort.

## Important caveats (read before trusting any number)

1. **`vaf_truth` is the expected panel-mean VAF, approximated as
 `0.5 x mixture_fraction_A`** (heterozygous-variant assumption; most HCC1395
 somatic SNVs are het). For precise quantification-bias, recompute per sample
 as the mean SEQC2-truth-set VAF over your panel variants times the fraction.
2. **Input mass is not constant** across samples. Different input = different
 unique-molecule depth, which affects LoD. Subset to one `input` value or treat
 input as a covariate.
3. **ILM2 is a standard ctDNA assay, NOT allele-specific enrichment.** Score with
 **enrichment factor = 1** (no de-biasing). This cohort validates consensus,
 interrogation, the empirical null, and the LoD/LoQ math, but CANNOT validate
 the 2Strands enrichment layer - that still needs in-house dilutions through the
 real probes.
4. **Cohort hygiene.** Split the 12 Bf blanks: fit the background/PoN on one half,
 evaluate FPR on the held-out half. `validate.py` does the disjoint split.

## How to fill it in (the real driver)

`bin/run_validation_chain.sh` automates the whole cohort end-to-end. Because the
box egress is throttled (~1-2 MB/s) and full runs are ~40 GB each, it downloads a
**downsampled** slice per run (`fetch_and_interrogate.sh -N 3000000`, i.e.
`fastq-dump -X`):

```bash
PATH=/data/alvin/envs/mrd/bin:$PATH bash bin/run_validation_chain.sh
# per run: fastq-dump -X 3M -> bwa-mem2 align -> interrogate.py
# then: build_background.py (from the Bf blank rows) -> background.tsv + pon.tsv
# validate.py run-real -> results/seqc2_validation.json
```

> **Downsampling caveat (important):** 3M spots against a genome-wide panel gives
> very sparse per-site coverage (~median depth 1). This proves the cfDNA pipeline
> **plumbing** on real reads end-to-end, but is **not** a real limit of detection.
> Full-depth runs, or a targeted-panel cfDNA cohort (coverage concentrated on
> panel sites), are needed for true LoD/LoQ.

To build the panel this cohort is interrogated against, see `bin/build_panel.sh`
(SEQC2 high-confidence HCC1395 truth INTERSECT confident regions; truth from
PRJNA489865). See [../../VALIDATION.md](../../VALIDATION.md) for the full runbook.

## Provenance

Built from NCBI E-utilities runinfo + BioSample attributes for PRJNA677999.
Mixture fractions are quoted verbatim from the BioSample `isolate` attributes.
Download each run with `prefetch <accession>`; no hard-coded URLs (SRA paths are
region/version specific).
