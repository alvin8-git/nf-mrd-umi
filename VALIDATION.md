# Validation Runbook

How to take this pipeline out of the silo: run it against real public data and
get real accuracy, specificity, and (eventually) LoD/LoQ numbers instead of
simulated ones. Every dataset and accession here was verified against NCBI/EGA
(see [VALIDATION_DATASETS.md](VALIDATION_DATASETS.md)); the per-component coverage
map and the full benchmark tables live in
[VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md).

This runbook now documents procedures that have actually been run, not just
planned ones:

- **Pipeline A (panel_design)** - run end-to-end on real SEQC2 HCC1395 tumor/normal
 WES and benchmarked against the SEQC2 gold truth set. **Done.**
- **Pipeline B (mrd_monitor)** - driven over the real SEQC2 ILM2 cfDNA titration
 cohort, downsampled. **In progress; proves plumbing, not LoD.**
- **Offline math + wiring** - the MRD statistics validated on simulated cohorts,
 and the Nextflow DAG validated with `-stub-run`. **Done, reproducible anytime.**

## What this validates (and what it does not)

Validates: the somatic-calling / panel-selection path (Pipeline A) against a
community gold standard; the cfDNA interrogation plumbing (Pipeline B) on real
reads end-to-end; and the empirical null, LoB/LoD/LoQ statistics on simulated
cohorts.

Does NOT validate yet: a real cfDNA limit-of-detection (the cohort is downsampled - 
see B), and the allele-specific enrichment layer. SEQC2 uses hybrid/amplicon
assays, not 2Strands probes, so score with enrichment = 1. The enrichment
de-biasing and per-probe calibration still need in-house dilutions through the
real probes.

---

## A. Pipeline A - real WES run + accuracy benchmark (DONE)

Run the panel-design pipeline on a matched tumor + normal WES pair and score the
somatic SNV calls against the SEQC2 truth set. We used SEQC2 **HCC1395** (tumor,
`SRR7890850`) + **HCC1395BL** (normal, `SRR7890851`), GRCh38.

### A.1 Download the WES pair

Resumable prefetch + fasterq-dump; drops the `.sra` after extracting fastq:

```bash
PATH=/data/alvin/envs/mrd/bin:$PATH \
 bash bin/fetch_wes.sh /data/alvin/wes SRR7890850 SRR7890851
# -> /data/alvin/wes/SRR7890850_{1,2}.fastq (tumor)
# /data/alvin/wes/SRR7890851_{1,2}.fastq (normal)
```

WES is far smaller than the 40 GB WGS runs and reuses the truth VCF already
downloaded for the panel.

### A.2 Stage the required resources

`bin/run_panel_design.sh` **fails loud** if a required input is missing, so nothing
silently no-ops. Already present on the box: the GRCh38 reference (`hg38.fa` +
`.fai` + `.dict`). You must supply two more, and may supply three best-practice
extras. Everything below is pulled from public buckets with `bin/fetch_url.py`
(a resumable urllib downloader - `curl`/`wget` are blocked in this env):

```bash
python3 bin/fetch_url.py <URL> <DEST> [<URL> <DEST> ...]
```

| env var | required | what | size |
|---|---|---|---|
| `SNP_VCF` | **yes** | chr-prefixed common-SNP VCF for FACETS pileup. We used GATK `1000G_phase1.snps.high_confidence.hg38.vcf.gz` (+ `.tbi`) | ~1.9 GB |
| `VEP_CACHE` | **yes** | ensembl-vep **116** GRCh38 offline cache directory | ~27.6 GB |
| `GERMLINE` | optional | Mutect2 germline resource `af-only-gnomad.hg38.vcf.gz` | large |
| `PON` | optional | Mutect2 panel-of-normals `1000g_pon.hg38.vcf.gz` | - |
| `INTERVALS` | optional | WES exome-target BED. **Without it Mutect2 scans the whole genome (~90 min).** | small |

Contigs must be `chr`-prefixed (the reference is UCSC-style). Egress is throttled
(~1-2 MB/s); the VEP cache is the long pole.

### A.3 Launch

```bash
SNP_VCF=/data/alvin/ref/GRCh38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
VEP_CACHE=/data/alvin/ref/vep_cache \
 bash bin/run_panel_design.sh
# add GERMLINE=... PON=... INTERVALS=... for the best-practice run (see script header)
```

Publishes `results/panel_design/HCC1395.panel.{bed,vcf,lock}`.

Smoke the DAG with no data first if you only want to check wiring:

```bash
nextflow run main.nf -profile docker,test -stub-run --workflow panel_design
```

### A.4 Result (this run)

- **50-site panel** selected.
- **FACETS purity 0.90 / ploidy 3.06** - matches the known near-pure aneuploid
 HCC1395 line (cross-check that the run is real, not stubbed; 33/50 selected sites
 are confirmed truth SNVs).
- **somatic-SNV F1 = 0.79** exome-restricted (precision 0.79 / recall 0.80).

### A.5 Benchmark method (reproducible)

1. Extract `FilterMutectCalls` PASS SNVs.
2. Reference-normalize: `bcftools norm -m-any -f <ref>`.
3. Match to the SEQC2 truth
 (`ref/seqc2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`) by
 `chrom:pos:ref>alt`.
4. Score in the **fair evaluation region** = SEQC2 HC regions intersect WES-callable
 (tumor depth >=10x). The callable footprint comes from
 `mosdepth --quantize 0:10:` on the tumor BAM (89.9 Mb).

The **genome-wide F1 (0.054) is an artifact** of benchmarking a WES (~1-2% of the
genome) against a whole-genome truth: ~38k truth SNVs sit in physically uncovered
space and become false FN. The **exome-restricted F1 (0.79) is the real number.**
0.79 is the expected ballpark for a single-caller, no-PoN, no-germline-resource
Mutect2 WES run; adding `GERMLINE`/`PON`/`INTERVALS` (A.2) is what closes the gap
to the ~0.85-0.95 published SEQC2 pipelines reach. The **full TP/FP/FN table is in
[VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md)** ("Pipeline A accuracy benchmark").

---

## B. Pipeline B - real cfDNA validation, downsampled (IN PROGRESS)

`bin/run_validation_chain.sh` drives the SEQC2 **ILM2** titration cohort headless:
download -> align -> interrogate per run, then build the background from the blank
rows, then score.

Cohort (manifest
[data/manifests/seqc2_ilm2_titration.manifest.tsv](data/manifests/seqc2_ilm2_titration.manifest.tsv),
**48 runs**): **12 blank** (0% tumor, `Bf`), **24 at 2%**, **12 at 10%**.

```
ILM2 titration (PRJNA677999)
 build_panel.sh (C) -> panel.vcf.gz truth SNVs INTERSECT HC regions
 fetch_and_interrogate -N 3e6 -> results/interrogate/<RUN>.tsv per-run counts
 build_background.py -> background.tsv, pon.tsv from the blank rows
 validate.py run-real -> results/seqc2_validation.json FPR/LoD/LoQ
```

### B.1 Run the chain

```bash
bash bin/run_validation_chain.sh
```

It waits out any in-flight single-run smoke driver (so workdirs don't collide),
runs the full downsampled batch (resumable - skips runs whose output already
exists), builds `background.tsv`/`pon.tsv` from the `blank` rows
(`build_background.py run --label blank`), and runs `validate.py run-real` ->
`results/seqc2_validation.json`.

The batch downsamples every run with `fetch_and_interrogate.sh -N 3000000`
(`fastq-dump -X` stops after the first 3M spots) because the box egress is
throttled (~1-2 MB/s) and full-depth runs are ~40 GB each.

### B.2 HONEST CAVEAT - this is plumbing, not LoD

Downsampling to **3M spots gives very sparse per-site coverage** (~median depth **1**
at the genome-wide panel sites). So this proves the cfDNA interrogation pipeline
runs end-to-end on **real reads** - download, align, dedup-flag, count - but it is
**NOT a real limit-of-detection.** A true LoD/LoQ needs either full-depth runs or a
**targeted-panel** cfDNA cohort, so coverage concentrates on panel sites instead of
spreading thin across the genome. (See the "Known gaps" section of
[VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md).)

### B.3 Read the result

`results/seqc2_validation.json` reports, per the held-out blanks and dilution
levels: observed FPR vs nominal alpha (0.05), LoD95, LoQ + per-level `rel_bias`, and
per-level hit-rate (sensitivity at 2% / 10%). Read these against the downsampling
caveat above - at median depth 1 they characterize the plumbing, not the assay.

---

## C. Build the SEQC2 panel (input to B)

The "panel" is the known HCC1395 mutations the cohort can see: SEQC2
high-confidence somatic SNVs, intersected with the confident regions (and,
optionally, an assay target BED). `bin/build_panel.sh` auto-downloads the verified
SEQC2 truth set + HC regions BED and intersects to `panel.vcf.gz` (+ `panel.bed`).

```bash
# genome-wide panel = truth SNVs INTERSECT HC regions
bin/build_panel.sh -o panel -r /data/alvin/ref/GRCh38/hg38.fa
# targeted assay: add -b ASSAY_TARGET.bed to also intersect the assay's regions
# add -C ucsc (or -C ensembl) if chrom naming differs from your reference
```

Truth set: SEQC2 Somatic Mutation WG (Fang et al. 2021), the same gold set
Pipeline A is benchmarked against. **Check the site count it prints** - if it warns
the panel is small (<8 sites), a hotspot assay BED may simply not cover enough
HCC1395 mutations for panel-integrated MRD. The output feeds
`interrogate.py --panel` (and the `PANEL=panel.vcf.gz` in `run_validation_chain.sh`).

---

## D. Offline validation - MRD math + wiring (DONE, reproducible)

These need no data or containers and should pass anytime.

**MRD statistics** (simulate blank + dilution cohorts, fit background on a training
fold, evaluate a disjoint held-out fold):

```bash
python3 bin/validate.py simulate-run --out /data/alvin/tmp/mrd_validation_report.json
# empirical null FPR ~0.008 (vs nominal 0.05); LoD95 VAF ~1e-4
```

**Nextflow wiring**, stub blocks only - independent of VEP cache / container
contents:

```bash
nextflow run main.nf -stub-run --workflow panel_design \
 --input assets/samplesheet_wes_hcc1395.csv ...
nextflow run main.nf -stub-run --workflow mrd_monitor \
 --input assets/samplesheet_cfdna.csv ...
```

The panel_design stub run passed (exit 0, all 14 tasks green) on the real HCC1395
samplesheet - see VALIDATION_COVERAGE.md "Pipeline A wiring validation".

> These are *simulated/stub* baselines. Real numbers worse than simulated is the
> expected, honest outcome: real data has artifacts the simulator does not model.
> That gap is the next round of work.

---

## Caveats (carried from the data, not bugs)

1. **Pipeline B is downsampled** (B.2): median depth ~1, so its FPR/LoD/LoQ
 characterize plumbing, not the assay. Full-depth or a targeted cfDNA cohort is
 needed for a real limit-of-detection.
2. **Pipeline A is single-caller** with no PoN / germline resource by default
 (F1 0.79). The `GERMLINE`/`PON`/`INTERVALS` knobs (A.2) are the path to ~0.9.
3. **Enrichment layer is not tested here** (no allele-specific probes in public
 data). Score with enrichment = 1.
4. **`vaf_truth` is approximate** (`0.5 x mixture fraction`, het assumption). Fine
 for LoD binning; recompute from the truth-set VAF distribution over your panel
 for precise quantification bias.
5. **No-UMI molecule counting** inflates depth (duplicates flagged, not collapsed).
 Use the UMI duplex path for true molecule-level counts on a UMI assay.
6. **Cohort hygiene** - never fit the background and test specificity on the same
 blanks; keep the folds disjoint.

## See also

- [README.md](README.md) - pipeline overview and the one-command self-tests
- [VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md) - per-component coverage map +
 the full Pipeline A benchmark and wiring-validation tables
- [VALIDATION_DATASETS.md](VALIDATION_DATASETS.md) - full dataset catalogue
- [data/manifests/README.md](data/manifests/README.md) - the cohorts and their caveats
- [TODO.md](TODO.md) - open validation items
- [Documentation.md](Documentation.md) - pipeline 101
