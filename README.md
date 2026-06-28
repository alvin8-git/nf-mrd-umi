# nf-mrd-umi

A tumor-informed **Minimal Residual Disease (MRD)** pipeline for an
allele-specific enrichment circulating-tumor-DNA (ctDNA) liquid-biopsy assay.
It detects cancer recurrence from a blood draw by tracking a patient's known
tumor mutations down to ultra-low frequencies (below 0.1% variant allele
fraction, toward 0.01%).

> **Status: research-stage, not for clinical use.** Both Nextflow pipelines are
> built (DSL2) and **Pipeline A is validated end-to-end on real SEQC2 data**
> (somatic-calling F1 = 0.79 exome-restricted; FACETS purity/ploidy match the
> HCC1395 truth). The MRD-call engine is self-tested and validated on a held-out
> simulated cohort; real-data MRD LoD is in progress. See [Status](#status).

New here? Read **[Documentation.md](Documentation.md)** first: it is a 101 on the
clinical problem, the vocabulary, the architecture, the tools, and every design
decision.

---

## What it is

Two Nextflow DSL2 pipelines (`--workflow panel_design | mrd_monitor`):

- **Pipeline A (panel design, once per patient):** matched tumor + buffy-coat
 whole-exome sequencing (WES) -> bwa-mem2 align -> Mutect2 somatic discovery ->
 FACETS purity/ploidy + PyClone-vi clonality -> CHIP/germline/probe filtering ->
 personalized panel of trackable mutations (`*.panel.{bed,vcf}`) + a patient-lock
 provenance token (`*.panel.lock`).
- **Pipeline B (MRD monitoring, every blood draw):** cfDNA + fgbio duplex UMI
 consensus -> verify the patient-lock -> interrogate the known panel sites ->
 panel-integrated, enrichment-de-biased MRD call (positive / negative /
 indeterminate, with tumor fraction, p-value, CI, and provenance stamp).

Standard tools (fgbio, bwa-mem2, GATK Mutect2, FACETS, PyClone-vi, VEP, Picard)
do the standard jobs. The MRD-specific logic - the part no off-the-shelf tool
provides - lives in the custom Python engine in `bin/`.

---

## Repository layout

```
nf-mrd-umi/
+-- README.md / Documentation.md / TODO.md # this / the 101 / status & gaps
+-- VALIDATION.md / VALIDATION_COVERAGE.md # runbook / component->dataset+results
+-- VALIDATION_DATASETS.md / CONTAINERS.md # public data / container strategy
+-- main.nf, nextflow.config, conf/ # DSL2 entry + profiles (docker/test/local)
+-- modules/local/ subworkflows/local/ workflows/ # 20 modules, 3 subworkflows, 2 workflows
+-- bin/ # the custom engine + glue + drivers
| +-- panel_select.py interrogate.py mrd_integrate.py # core MRD engine
| +-- validate.py sample_id.py build_background.py # QC: harness / identity / null
| +-- pyclone_prep.py normal_evidence.py # Pipeline A glue
| +-- run_panel_design.sh run_validation_chain.sh ... # real-run drivers
+-- assets/example_panel_HCC1395/ # a real validated Pipeline A panel
+-- .github/workflows/ci.yml # selftests + both-pipeline stub-runs
```

---

## Quickstart

**Every self-test (8 engine scripts), one command:**

```bash
bash bin/run_selftests.sh # panel_select, interrogate, mrd_integrate, validate,
 # build_background, pyclone_prep, normal_evidence, sample_id
```

**Validate the MRD math** (simulates blank + dilution cohorts, fits background on a
training fold, evaluates a disjoint held-out fold, reports specificity/sensitivity):

```bash
python3 bin/validate.py simulate-run --out /tmp/mrd_validation_report.json
# empirical null FPR ~0.008 (vs nominal 0.05); LoD95 VAF ~1e-4
```

**Validate the Nextflow wiring offline** (no data, no containers - stub blocks only):

```bash
nextflow run main.nf -stub-run --workflow panel_design --input assets/samplesheet_wes_hcc1395.csv ...
nextflow run main.nf -stub-run --workflow mrd_monitor --input assets/samplesheet_cfdna.csv ...
```

**Run Pipeline A for real** (needs reference + resources; see `bin/run_panel_design.sh`
header and [VALIDATION.md](VALIDATION.md)):

```bash
bash bin/run_panel_design.sh # tumor/normal WES -> results/panel_design/*.panel.{bed,vcf,lock}
```

---

## The custom engine

| Script | Pipeline / role | Commands |
|---|---|---|
| `panel_select.py` | A - personalized panel selection (CHIP/germline/probe/CCF) | `run`, `selftest` |
| `interrogate.py` | B - duplex-aware unique-molecule counting at panel sites | `run`, `selftest` |
| `mrd_integrate.py` | B - panel-integrated MRD call (empirical null + enrichment de-bias + provenance) | `run`, `selftest` |
| `validate.py` | QC - offline harness (FPR, LoB, LoD95, LoQ, enrichment recovery) | `simulate-run`, `run-real`, `selftest` |
| `sample_id.py` | QC - SNP-fingerprint concordance + patient-lock token (sample-swap protection) | `fingerprint`, `concordance`, `provenance`, `verify-lock`, `selftest` |
| `build_background.py` | QC - healthy-donor blanks -> beta-binomial background + PoN | `run`, `selftest` |
| `pyclone_prep.py`, `normal_evidence.py` | A glue - PyClone I/O + buffy-coat CHIP evidence | `build`/`to-ccf`, `pileup` |

---

## Status

| Component | State |
|---|---|
| Custom MRD engine (`bin/`, 8 self-tested scripts) | Implemented, self-tested |
| Empirical covariance-preserving null | Implemented, validated (FPR 0.28 -> 0.07) |
| Odds-space enrichment de-biasing | Implemented, validated (recovery 0.43x -> 1.0x) |
| **Pipeline A (panel_design) DSL2** | **Built + validated on real SEQC2 HCC1395 WES** |
| **Pipeline B (mrd_monitor) DSL2** | **Built; engine validated; real-data LoD in progress** |
| **Sample identity / patient-lock / provenance** | **Implemented (`sample_id.py`) and wired into both pipelines** |
| CI (selftests + stub-runs), PostToolUse hook, container pins | Implemented |
| Run controls, full-depth cfDNA LoD, multi-caller consensus | Not yet (see TODO.md) |
| SOPs / regulatory docs | Not started |

**Real-data result (Pipeline A, SEQC2 HCC1395 / HCC1395BL WES):** 50-site panel,
FACETS purity 0.90 / ploidy 3.06 (matches truth), somatic-SNV **F1 = 0.79**
(precision 0.79 / recall 0.80) in the exome-callable HC region. A re-run with a
Mutect2 panel-of-normals + gnomAD germline resource (toward ~0.9) is underway.
Full method and numbers: [VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md). The
validated panel is committed at `assets/example_panel_HCC1395/`.

See [TODO.md](TODO.md) for remaining gaps and [VALIDATION.md](VALIDATION.md) for
the real-world validation runbook.

---

## Disclaimer

This software is for research and development only. It is not a validated
diagnostic and must not be used for clinical decision-making in its current form.
