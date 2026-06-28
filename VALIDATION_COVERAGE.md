# Validation coverage — which public data exercises which pipeline component

The SEQC2 ILM2 titration (PRJNA677999) currently running validates only **Pipeline B's
no-UMI interrogation + VAF-vs-dilution** path. Below is the map of every *other*
component to a public dataset that would validate it. All accessions must be
re-verified to resolve before download (egress is throttled to ~1-2 MB/s here).

| # | Component (unvalidated) | Dataset / accession | What it proves | Caveat |
|---|---|---|---|---|
| 0 | Pipeline B no-UMI interrogate + VAF | **SEQC2 ILM2 titration** PRJNA677999 (running) | download->align->interrogate->VAF vs known dilution (0/2/10%) | downsampled = plumbing only, not LoD (median depth 1) |
| 1 | **Pipeline A: somatic calling + filtering** (Mutect2 -> FilterMutectCalls -> VEP) | **HCC1395 tumor + HCC1395BL normal WES/WGS**, PRJNA677999 | precision/recall of called somatic SNV/indel vs the SEQC2 **high-confidence truth VCF** (the same gold set we build the panel from) | highest-value add: closes Pipeline A end-to-end against a community gold standard |
| 2 | **Sample identity / patient-lock** (`bin/sample_id.py`) | HCC1395 + HCC1395BL (same donor) **+ NA12878 / HG001** (GIAB, unrelated) | concordance(HCC1395, HCC1395BL) -> SAME; concordance(HCC1395, NA12878) -> DIFFERENT — real swap detection | needs a common-SNP sites file (GATK fingerprint HaplotypeMap, or a gnomAD common-SNP subset) |
| 3 | **UMI duplex consensus** (fgbio subworkflow) | **nf-core/fastquorum tiny `SRR6109255`** (~2.7 MB; full run 276 MB), read structure **`10M1S+T 10M1S+T`** — CRISPR-DS duplex of TP53, MiSeq | FastqToBam->GroupReadsByUmi->CallDuplexConsensus->Filter; duplex_support > 0 | read structure IS documented (nf-core test-datasets `fastquorum` branch). For the cfDNA CorrectUmis path: IDT xGen sheet, `8M+T 8M+T` |
| 4 | **FACETS CN + PyClone CCF** (clonal structure + `pyclone_prep`) | HCC1395 tumor/normal, PRJNA677999 | FACETS purity/ploidy + PyClone clustering run to completion; plausible truncal CCF feeds panel_select | soft validation — HCC1395 CN is characterized but there is no exact per-variant CCF truth |
| 5 | **Empirical healthy-donor cfDNA null** (background) | in-cohort: ILM2 **Bf blanks** (0% tumor, already used); external: a true healthy-donor cfDNA cohort | covariance-preserving empirical null calibrated on real blanks | most external healthy-cfDNA cohorts are controlled-access (EGA/dbGaP) |
| 6 | **Allele-specific enrichment calibration** | **none public** — proprietary to the 2Strands assay | odds-space de-bias multiplier from a dilution series | can only be validated on the real assay's dilution series; standard libs have enrich=1.0 |

## Highest-leverage next download (recommendation)

**#1 — SEQC2 HCC1395 / HCC1395BL WES (PRJNA677999).** One tumor/normal pair unlocks
THREE validations at once:
- Pipeline A somatic calling vs the SEQC2 gold truth (precision/recall).
- Sample identity (tumor vs its own matched normal -> SAME).
- FACETS + PyClone on real clonal structure.

It is the same BioProject already in use, WES is far smaller than the 40 GB WGS runs,
and it reuses the truth VCF already downloaded for the panel. Add one unrelated GIAB
sample (NA12878) to complete the identity DIFFERENT case.

## Known gaps that no public data closes
- Enrichment calibration (#6) — assay-proprietary.
- A large matched healthy-donor cfDNA null (#5 external) — controlled-access.
- Full-depth cfDNA LoD/LoQ — egress-bound on this box; needs the full-depth runs or a
  targeted (not WGS) cfDNA cohort so coverage concentrates on panel sites.

## Pipeline A wiring validation (2026-06-28)

`nextflow run main.nf -stub-run --workflow panel_design` with the real HCC1395 WES
samplesheet + resource paths **passed (exit 0)** — all 14 tasks green:

```
ALIGN_WES (BWAMEM2_MEM T+N -> SAMTOOLS_SORT -> PICARD_MARKDUPLICATES)
GATK4_MUTECT2 -> GATK4_FILTERMUTECTCALLS
FACETS | VEP
PYCLONE_PREP -> PYCLONEVI
NORMAL_EVIDENCE
PANEL_SELECT -> HCC1395.panel.bed + HCC1395.panel.vcf
```

Confirms tumor/normal pairing by patient and the three-way joins
(`somatic x FACETS.cnv x purity`, `VEP x ccf x normal_evidence`) resolve. `-stub-run`
exercises wiring only (stub touch/echo blocks), so it is independent of VEP cache
contents. Launch the real run via `bin/run_panel_design.sh` once the VEP cache finishes
and the box egress frees up (first run also pulls the gatk4/picard/vep/facets/pyclone-vi
biocontainers pinned in `conf/docker.config`).
