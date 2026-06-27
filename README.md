# nf-mrd-umi

A tumor-informed **Minimal Residual Disease (MRD)** pipeline for an
allele-specific enrichment circulating-tumor-DNA (ctDNA) liquid-biopsy assay.
It detects cancer recurrence from a blood draw by tracking a patient's known
tumor mutations down to ultra-low frequencies (below 0.1% variant allele
fraction, toward 0.01%).

> **Status: research-stage, not for clinical use.** The custom statistical
> engine (`bin/`) is implemented and self-tested. The Nextflow orchestration is
> designed but not yet built. See [Status](#status) and [TODO.md](TODO.md).

New here? Read **[Documentation.md](Documentation.md)** first: it is a 101 on the
clinical problem, the vocabulary, the architecture, the tools, and every design
decision.

---

## What it is

Two pipelines:

- **Pipeline A (panel design, once per patient):** matched tumor + buffy-coat
  whole-exome sequencing (WES) -> somatic discovery -> purity/ploidy/clonality ->
  CHIP/germline filtering -> personalized panel of trackable mutations.
- **Pipeline B (MRD monitoring, every blood draw):** cfDNA + duplex UMI consensus
  -> interrogate the known panel sites -> panel-integrated, enrichment-de-biased
  MRD call (positive / negative / indeterminate, with tumor fraction, p-value,
  and confidence interval).

The standard genomics tools (fgbio, bwa-mem2, GATK Mutect2, PURPLE/FACETS,
PyClone-vi, VEP) do the standard jobs. The MRD-specific logic, the part no
off-the-shelf tool provides, lives in four custom Python scripts.

---

## Repository layout

```
nf-mrd-umi/
+-- README.md            # this file
+-- Documentation.md     # the 101: problem, terms, architecture, design decisions
+-- TODO.md              # roadmap and known gaps
+-- CLAUDE.md            # agent/tooling routing rules
+-- bin/                 # the custom MRD engine (implemented, self-tested)
    +-- panel_select.py    # Pipeline A: personalized panel selection
    +-- interrogate.py     # Pipeline B: molecule counting at panel sites
    +-- mrd_integrate.py   # Pipeline B: the MRD detection engine
    +-- validate.py        # offline validation harness (LoB/LoD/LoQ, FPR)
```

The Nextflow DSL2 modules and workflows are specified in the design docs but are
not yet committed as `.nf` files.

---

## Requirements

The custom engine (everything in `bin/`):

- Python 3.10+
- `numpy`, `scipy` (statistics), `pysam` (only for the BAM/VCF I/O paths;
  the `selftest` subcommands run without it)

The full pipeline additionally needs (once the Nextflow layer is built):
Nextflow, fgbio, bwa-mem2, GATK, CNVkit, PURPLE/FACETS/TITAN, PyClone-vi,
Ensembl VEP, Picard, and the GRCh38 no-alt + decoy reference.

---

## Quickstart

Everything in `bin/` is runnable today on synthetic data, no sequencing required.

**Run every self-test:**

```bash
python3 bin/panel_select.py selftest
python3 bin/interrogate.py selftest
python3 bin/mrd_integrate.py selftest
python3 bin/validate.py selftest
```

**Run the validation harness** (simulates blank + dilution cohorts, fits the
background on a training fold, evaluates a disjoint held-out fold, and reports
specificity and sensitivity):

```bash
python3 bin/validate.py simulate-run --out /tmp/mrd_validation_report.json
```

Example output:

```
  null model        : empirical_donor_resampling
  observed FPR      : 0.008   (OK)            # vs nominal alpha 0.05
  LoB (tumor frac)  : 7.78e-06
  LoD95 VAF         : 0.0001
  LoQ VAF           : 0.0001
  enrichment recov. : 1.004x of truth
```

**Compare the empirical null against the old independent null** on identical
data (demonstrates why the empirical null matters):

```bash
python3 bin/validate.py simulate-run --null independent --out /tmp/indep.json
python3 bin/validate.py simulate-run --null empirical   --out /tmp/emp.json
# independent FPR ~0.28 (inflated) vs empirical ~0.01 (controlled)
```

---

## The four scripts

| Script | Pipeline / stage | Command | Does |
|---|---|---|---|
| `panel_select.py` | A - panel design | `run`, `selftest` | Annotated somatic VCF -> personalized panel (BED + VCF), with CHIP/germline/probe filters and clonality ranking |
| `interrogate.py` | B - stage B5 | `run`, `selftest` | Consensus BAM + panel -> per-site unique-molecule counts (duplex-aware) |
| `mrd_integrate.py` | B - stage B6 | `run`, `selftest` | Per-site counts + background -> MRD call, tumor fraction, p-value, CI |
| `validate.py` | offline QC | `simulate-run`, `run-real`, `selftest` | Measures false-positive rate, LoB, LoD95, LoQ, enrichment recovery |

Each `run`/`simulate-run` command prints its own `--help`. Example end-to-end
(Pipeline B scoring of one timepoint):

```bash
python3 bin/interrogate.py run \
    --bam consensus.bam --panel panel.vcf --out site_counts.tsv

python3 bin/mrd_integrate.py run \
    --site-counts site_counts.tsv \
    --background healthy_donor_background.tsv \
    --pon healthy_donor_pon.tsv \
    --patient-id PT001 --timepoint T1 \
    --out PT001_T1.mrd.json
```

---

## Status

| Component | State |
|---|---|
| `bin/` custom engine (4 scripts) | Implemented, self-tested |
| Empirical covariance-preserving null | Implemented, validated (FPR 0.28 -> 0.07) |
| Odds-space enrichment de-biasing | Implemented, validated (recovery 0.43x -> 1.0x) |
| Validation harness | Implemented (`simulate-run`, `run-real`, `selftest`) |
| Nextflow DSL2 modules / workflows | Designed, not yet built |
| Sample identity, run controls, provenance | Not yet implemented (required for clinical use) |
| SOPs / regulatory docs | Not started |

**Self-tests prove the code implements its model, not that the model matches
biology.** Real validity requires `validate.py run-real` on wet-lab dilution
series and healthy-donor cohorts.

See [TODO.md](TODO.md) for the full roadmap (P0 clinical scaffolding, P1
specificity hardening, P2 sensitivity and scope), and
[VALIDATION.md](VALIDATION.md) for the real-world validation runbook (SEQC2
public data -> real FPR/LoD/LoQ).

---

## Disclaimer

This software is for research and development only. It is not a validated
diagnostic and must not be used for clinical decision-making in its current form.
