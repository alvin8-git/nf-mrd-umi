# Validation coverage  -  which public data exercises which pipeline component

The SEQC2 ILM2 titration (PRJNA677999) currently running validates only **Pipeline B's
no-UMI interrogation + VAF-vs-dilution** path. Below is the map of every *other*
component to a public dataset that would validate it. All accessions must be
re-verified to resolve before download (egress is throttled to ~1-2 MB/s here).

| # | Component (unvalidated) | Dataset / accession | What it proves | Caveat |
|---|---|---|---|---|
| 0 | Pipeline B no-UMI interrogate + VAF | **SEQC2 ILM2 titration** PRJNA677999 (running) | download->align->interrogate->VAF vs known dilution (0/2/10%) | downsampled = plumbing only, not LoD (median depth 1) |
| 1 | **Pipeline A: somatic calling + filtering** (Mutect2 -> FilterMutectCalls -> VEP) | **HCC1395 tumor + HCC1395BL normal WES/WGS**, PRJNA677999 | precision/recall of called somatic SNV/indel vs the SEQC2 **high-confidence truth VCF** (the same gold set we build the panel from) | highest-value add: closes Pipeline A end-to-end against a community gold standard |
| 2 | **Sample identity / patient-lock** (`bin/sample_id.py`) | HCC1395 + HCC1395BL (same donor) **+ NA12878 / HG001** (GIAB, unrelated) | concordance(HCC1395, HCC1395BL) -> SAME; concordance(HCC1395, NA12878) -> DIFFERENT  -  real swap detection | needs a common-SNP sites file (GATK fingerprint HaplotypeMap, or a gnomAD common-SNP subset) |
| 3 | **UMI duplex consensus** (fgbio subworkflow) | **nf-core/fastquorum tiny `SRR6109255`** (~2.7 MB; full run 276 MB), read structure **`10M1S+T 10M1S+T`**  -  CRISPR-DS duplex of TP53, MiSeq | FastqToBam->GroupReadsByUmi->CallDuplexConsensus->Filter; duplex_support > 0 | read structure IS documented (nf-core test-datasets `fastquorum` branch). For the cfDNA CorrectUmis path: IDT xGen sheet, `8M+T 8M+T` |
| 4 | **FACETS CN + PyClone CCF** (clonal structure + `pyclone_prep`) | HCC1395 tumor/normal, PRJNA677999 | FACETS purity/ploidy + PyClone clustering run to completion; plausible truncal CCF feeds panel_select | soft validation  -  HCC1395 CN is characterized but there is no exact per-variant CCF truth |
| 5 | **Empirical healthy-donor cfDNA null** (background) | in-cohort: ILM2 **Bf blanks** (0% tumor, already used); external: a true healthy-donor cfDNA cohort | covariance-preserving empirical null calibrated on real blanks | most external healthy-cfDNA cohorts are controlled-access (EGA/dbGaP) |
| 6 | **Allele-specific enrichment calibration** | **none public**  -  proprietary to the 2Strands assay | odds-space de-bias multiplier from a dilution series | can only be validated on the real assay's dilution series; standard libs have enrich=1.0 |

## Highest-leverage next download (recommendation)

**#1  -  SEQC2 HCC1395 / HCC1395BL WES (PRJNA677999).** One tumor/normal pair unlocks
THREE validations at once:
- Pipeline A somatic calling vs the SEQC2 gold truth (precision/recall).
- Sample identity (tumor vs its own matched normal -> SAME).
- FACETS + PyClone on real clonal structure.

It is the same BioProject already in use, WES is far smaller than the 40 GB WGS runs,
and it reuses the truth VCF already downloaded for the panel. Add one unrelated GIAB
sample (NA12878) to complete the identity DIFFERENT case.

## Known gaps that no public data closes
- Enrichment calibration (#6)  -  assay-proprietary.
- A large matched healthy-donor cfDNA null (#5 external)  -  controlled-access.
- Full-depth cfDNA LoD/LoQ  -  egress-bound on this box; needs the full-depth runs or a
  targeted (not WGS) cfDNA cohort so coverage concentrates on panel sites.

## Pipeline A wiring validation (2026-06-28)

`nextflow run main.nf -stub-run --workflow panel_design` with the real HCC1395 WES
samplesheet + resource paths **passed (exit 0)**  -  all 14 tasks green:

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

## Pipeline A accuracy benchmark vs SEQC2 truth (2026-06-28)

Ran Pipeline A end-to-end on real SEQC2 WES (tumor HCC1395 SRR7890850 / normal
HCC1395BL SRR7890851, GRCh38) and benchmarked the somatic SNV calls against the
SEQC2 high-confidence truth set (`high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`).

**Method.** PASS SNVs from FilterMutectCalls, reference-normalized (`bcftools norm
-m-any -f`), matched to truth by `chrom:pos:ref>alt`. Two evaluation regions:
- *genome-wide* = SEQC2 HC regions (the truth's native domain).
- *exome-callable* = HC INT tumour depth >=10x, where the WES could physically call.
  Callable footprint from `mosdepth --quantize 0:10:` on the tumour BAM (89.9 Mb).

| metric | genome-wide (artifact) | exome-callable (fair) |
|---|---|---|
| TP | 1103 | 1032 |
| FP | 334  | 281  |
| FN | 38344 | 265 |
| Precision | 0.77 | **0.79** |
| Recall | 0.028 | **0.80** |
| **F1** | 0.054 | **0.79** |

The genome-wide F1 (0.054) is an **artifact**: WES covers ~1-2% of the genome but the
truth is genome-wide, so ~38k truth SNVs are physically uncovered -> false FN. Once
the truth is restricted to the WES-callable space, recall is 0.80 and **F1 = 0.79**.

**Cross-checks that the run is real, not stubbed:** FACETS purity **0.90** / ploidy
**3.06** match the known near-pure aneuploid HCC1395 line; 33/50 selected panel sites
are confirmed truth SNVs.

**Interpretation.** 0.79 is the expected ballpark for a **single-caller, no-PoN,
no-germline-resource** Mutect2 WES run (all optional and unset here). Published SEQC2
WES pipelines reach ~0.85-0.95 using a Mutect2 panel-of-normals + gnomAD germline
resource + multi-caller consensus. The FP=281 (germline/artifact leakage) and FN=265
(subclonal/low-VAF/filtered) are exactly what those additions address; the knobs are
already wired in `bin/run_panel_design.sh` (`GERMLINE=`, `PON=`, `INTERVALS=`).

### Re-run with Mutect2 PoN + gnomAD (2026-06-28)

Re-ran Pipeline A with the GATK `1000g_pon.hg38` panel-of-normals + `af-only-gnomad.hg38`
germline resource (same exome-callable eval region). The germline/artifact filtering
did what it should:

| metric (exome-callable HC region) | baseline (none) | **PoN + gnomAD** |
|---|---|---|
| PASS SNVs (total) | 2877 | 1736 |
| TP | 1032 | 1024 |
| FP | 281 | **221** (-60, ~21%) |
| FN | 265 | 273 |
| Precision | 0.786 | **0.822** |
| Recall | 0.796 | 0.790 |
| **F1** | 0.791 | **0.806** |

The bigger win is the **panel** (the deliverable, not the full call set): **50/50**
selected sites are now confirmed SEQC2 truth somatic SNVs, up from **33/50** without
PoN+gnomAD - the 17 germline/artifact contaminants are gone, so every tracked site is
a bona-fide tumor mutation. The committed example (`assets/example_panel_HCC1395/`) is
this cleaner panel.

**Decision (matched PoN):** not pursued. The generic 1000G PoN already captured the
achievable FP reduction; a center-matched HCC1395BL PoN would cost ~228 GB download +
12 Mutect2 runs, is cross-kit (only 1 LL normal among 6 centers), and the remaining 221
FPs are largely SEQC2-truth-conservative. The real levers toward ~0.9 F1 are multi-caller
consensus (add Strelka2) and the FN/subclonal side, not a better PoN.
