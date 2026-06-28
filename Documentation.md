# nf-mrd-umi: Pipeline 101

A from-scratch background to this pipeline: the clinical problem, the vocabulary,
the architecture, the tools, and the reasons behind each design decision. If you
are new to the project (or to circulating-tumor-DNA assays), start here. For how
to actually run things, see [README.md](README.md). For what is left to build,
see [TODO.md](TODO.md).

---

## 1. What this pipeline does, in one paragraph

It detects whether a cancer patient's tumor is coming back, from a tube of blood,
before a scan could ever see it. It does this by first reading the patient's
tumor (from a tissue biopsy) to learn that tumor's specific DNA mutations, then
watching for those exact mutations in the blood over time. A positive signal at
extremely low levels (below 1 mutant molecule in 1,000, sometimes below 1 in
10,000) is an early warning of residual or returning disease. This is called
**Minimal Residual Disease (MRD)** monitoring, and this is a **tumor-informed**
MRD pipeline built around an **allele-specific enrichment** assay.

---

## 2. The clinical problem (and why it is hard)

When a tumor grows, it sheds fragments of its DNA into the bloodstream. These
fragments are called **circulating tumor DNA (ctDNA)**, and they float in the
blood alongside the normal cell-free DNA (**cfDNA**) shed by healthy tissue.
After a patient is treated, a small amount of tumor may survive. If we could
detect its ctDNA, we could catch recurrence months earlier than imaging.

The difficulty is the signal-to-noise ratio:

- After successful treatment, ctDNA may be **less than 0.1%** of the cfDNA at a
 given mutation site, and sometimes **below 0.01%**. That is one mutant molecule
 hiding among ten thousand normal ones.
- Sequencing machines make mistakes at roughly the same rate (around 0.1% per
 base). So the true signal is buried at or below the noise floor of the
 instrument.
- Biology adds more noise: aging blood cells accumulate their own mutations
 (see **CHIP**, Section 4) that look exactly like tumor signal if you are not
 careful.

Beating this needs three things working together: a way to **remove sequencing
errors** (UMIs and consensus calling), a way to **concentrate the true signal**
(allele-specific enrichment), and a way to **integrate weak evidence across many
mutations** so that no single site has to carry the call (panel integration).
Those three ideas are the spine of this pipeline.

---

## 3. How the assay is validated: public ctDNA data vs the 2Strands probes

A pipeline that claims sub-0.1% detection has to be *proven*, not asserted. Two
kinds of data prove different parts of it, and it helps to see early why neither
alone is enough.

**The shape of a validation cohort.** Whatever the assay, you validate an MRD
pipeline with the same cohort shape:

- **Blanks** (tumor-free cfDNA): no tumor signal, so anything the pipeline calls
  in them is noise. The blanks *are* the background/null (`build_background.py`),
  and held-out blanks measure specificity / Limit of Blank.
- **A dilution series** (tumor titrated into normal at known fractions): each
  level is a sample whose true VAF you already know, so you can ask "did we detect
  it, and was the reported tumor fraction right?" -- that gives Limit of
  Detection, Limit of Quantitation, and quantification bias.

**What the public SEQC2 ILM2 cohort gives us.** The 48-run ILM2 titration
(HCC1395 tumor into HCC1395BL normal: 12 blanks at 0%, 24 at 2%, 12 at 10%) is
exactly that shape, on a *standard* Illumina ctDNA assay. It validates consensus,
interrogation, the empirical null, and the LoB/LoD/LoQ math on real reads. Because
HCC1395 is also the cell line Pipeline A discovers mutations from, the two halves
meet on the same "patient": Pipeline A finds the mutations, the ILM2 cfDNA
dilutions are where Pipeline B tracks them. (See
[data/manifests/README.md](data/manifests/README.md) and
[VALIDATION_DATASETS.md](VALIDATION_DATASETS.md).)

**What it CANNOT give us: the enrichment layer.** ILM2 is a standard assay, so its
panel sites read at roughly their true VAF (enrichment = 1). The 2Strands assay is
different: allele-specific probes *amplify the mutant allele*, so the observed VAF
is pushed up. A single real 2Strands run is, in role, like one ILM2 dilution
sample at a known VAF -- but its observed VAF is inflated, and recovering the true
VAF is the whole job.

**Recovering the true VAF when the enrichment factor is unknown.** For one unknown
sample you cannot: `observed_odds = enrichment x true_odds` is one equation with
two unknowns. So enrichment is not estimated per patient; it is *calibrated once*
as a property of the assay, then applied:

1. **Calibrate** on a 2Strands dilution series with known true VAFs. With true and
   observed both in hand, enrichment is the only unknown, and
   `fit_enrichment_from_dilution()` solves for it (panel-average, or per-probe if
   there is enough signal per site).
2. **Apply** to a patient sample: enrichment is now a known constant, so invert it
   -- `true_odds = observed_odds / enrichment` -- and read off the true VAF
   (`estimate_vaf(..., apply_enrichment=True)`).

It is a weighing scale: calibrate with known weights once, then weigh unknowns.

**Why odds-space.** Enrichment multiplies *odds* by a constant; in VAF-space the
same effect saturates (a 10x boost cannot take 30% to 300%). Odds-space keeps the
relationship linear and invertible. This is not cosmetic: the VAF-space fit
recovered only 0.43x of the true enrichment, the odds-space fit ~1.0x (see
[Section 7.3](#73-allele-specific-enrichment-and-odds-space-de-biasing) and the
`validate.py` numbers).

**The two things still required.** Calibrate-once-apply assumes enrichment is
*stable* (same multiplier in the patient run as at calibration) and acts as a
constant odds-multiplier. Both are checkable but unproven here:

- **Drift** is caught by a **per-run enrichment control** (a known-VAF spike-in
  every run) -- a run-level control on the roadmap ([TODO.md](TODO.md) P0).
- **Model correctness** -- that real probes behave as a clean odds-multiplier --
  can only be confirmed on a real 2Strands dilution series; the simulator cannot
  prove it (the open question in [TODO.md](TODO.md)). That low-VAF series (ideally
  down to 0.01-0.1%, where enrichment earns its keep) is the one piece no public
  dataset can substitute.


**In-silico stand-in (until the real series arrives).** `validate.py
simulate-from-real` builds exactly this cohort semi-synthetically: real ILM2
blanks supply the per-site depth/error structure, a synthetic forward enrichment
plus dilution plants healthy / tumor-ladder / a held-out dummy-patient sample, and
the real calibrate -> de-bias -> call flow runs on top. Its
`--enrichment-model constant|saturating` knob deliberately makes the forward model
*differ* from the de-biaser's constant-odds assumption, so it measures how
gracefully recovery degrades under model misspecification (enrichment recovery
~1.0x well-specified vs ~0.67x saturating; the dummy 0.05% patient is recovered to
within a few percent in the well-specified case). It validates the integration on
realistic structure -- it does **not** replace the real 2Strands dilution series.

### 3.1 A concrete walk-through (the ILM2 cohort through Pipeline B)

The pipeline is easiest to understand by following two real ILM2 runs through it
-- one blank, one tumor-spiked. `bin/run_validation_chain.sh` does all of this;
here is what each step produces.

1. **Fetch + align.** For each run accession (a 0% blank, and the 10% sample
   `SRR13385630`): `fastq-dump -X 3000000` pulls a downsampled slice and
   `bwa-mem2` aligns it to GRCh38 -> a sorted BAM.
2. **Interrogate the panel.** `interrogate.py run --bam <run>.bam --panel
   panel.vcf.gz` counts unique molecules at each of the 39,447 HCC1395 panel
   sites -> `<run>.site_counts.tsv` (per site: depth, alt-count, duplex support).
   In the blank, alt-counts at the tumor sites are ~0 (just error); in the 10%
   sample they sit near the expected VAF.
3. **Build the background from the blanks.** `build_background.py` fits a per-site
   beta-binomial error model from the 12 Bf blank tables -> `background.tsv`
   (+ `pon.tsv`, the donor-by-site matrix for the empirical null). This is the
   "what does noise look like here" model.
4. **Score each sample.** `mrd_integrate.py run --site-counts <run>.tsv
   --background background.tsv --pon pon.tsv` pools the panel sites into one
   tumor-fraction estimate, a Monte-Carlo p-value against the empirical null, and
   a bootstrap CI -> `<run>.mrd.json` (call: positive / negative / indeterminate).
5. **Read the cohort out.** `validate.py run-real` fits on a held-out split of
   blanks and reports the false-positive rate on the other blanks (specificity)
   and the hit-rate at each VAF level (2%, 10%).

Expected shape: blanks call negative (a few false positives bound specificity),
10% samples call positive comfortably, 2% samples sit near the limit of detection.
**Caveat:** at the downsampled depth used here (median ~1 molecule/site) those
numbers only demonstrate the flow; a real specificity/LoD needs full-depth runs
(Section 9, [VALIDATION.md](VALIDATION.md)).

## 4. The two-pipeline architecture

The work splits cleanly into two pipelines that run on different schedules and
different inputs. Keeping them separate is a deliberate decision (Section 6.1).

```
PIPELINE A - Panel Design (run ONCE per patient)
-------------------------------------------------
tumor biopsy WES -+
 +- align - somatic discovery - purity/ploidy/CN - clonality
buffy-coat WES --+ (Mutect2) (PURPLE/FACETS/CNVkit) (PyClone-vi)
 |
 annotate + germline/CHIP filter <-------+
 (VEP, gnomAD/dbSNP, buffy subtraction)
 |
 probe-design feasibility <--------------+
 (adapt 2Strands probe workflow)
 |
 PERSONALIZED PANEL ---------------------> (to Pipeline B)
 (BED + VCF + per-probe enrichment factors)

PIPELINE B - MRD Monitoring (run at EVERY blood draw / timepoint)
-------------------------------------------------
cfDNA blood - UMI extract - align - group by UMI - duplex consensus - realign
 (fgbio) (bwa-mem2) (fgbio)
 |
 interrogate panel sites <---------------+
 (count unique molecules supporting each known mutation)
 |
 panel-integrated call <------------------+
 (empirical null + enrichment de-bias)
 |
 MRD POSITIVE / NEGATIVE / INDETERMINATE
 + tumor fraction + p-value + confidence interval
```

**Pipeline A** answers "what should we watch for in this patient?" It reads the
tumor and the patient's healthy blood cells, finds the tumor's mutations, decides
which ones are reliable trackers, and designs the assay. It runs once.

**Pipeline B** answers "is the tumor here right now?" It takes a blood draw,
cleans up the sequencing errors, looks only at the mutations Pipeline A chose, and
produces a single yes/no/unsure call. It runs every time the patient gives blood.

---

## 5. Vocabulary (the terms you will keep meeting)

| Term | Plain meaning |
|---|---|
| **MRD** (Minimal Residual Disease) | Tiny amounts of cancer left after treatment, below imaging's reach. |
| **ctDNA** | Circulating tumor DNA: tumor-derived DNA fragments in the blood. |
| **cfDNA** | Cell-free DNA: all the free-floating DNA in blood, mostly from healthy cells. ctDNA is the tumor subset of cfDNA. |
| **Tumor-informed** | The assay is personalized: we already know the patient's mutations from their tumor, so we interrogate known sites instead of searching blind. |
| **VAF** (Variant Allele Fraction) | Fraction of molecules at a site carrying the mutant base. 0.1% VAF = 1 mutant in 1,000. |
| **Tumor fraction** | The overall fraction of cfDNA that is tumor-derived. For a heterozygous mutation, roughly 2x the VAF. |
| **UMI** (Unique Molecular Identifier) | A short random barcode attached to each original DNA molecule before PCR copying, so we can tell true molecules apart from PCR duplicates. |
| **Consensus read** | One error-corrected read rebuilt from all the PCR copies that share a UMI. One consensus read = one original molecule. |
| **Duplex** | Consensus that requires BOTH strands of the original double-stranded molecule to agree. The gold standard for error suppression. Not the same as paired-end (Section 6.2). |
| **Somatic vs germline** | Somatic = mutations the tumor acquired (what we want). Germline = inherited variants present in every cell (noise to remove). |
| **CHIP** (Clonal Hematopoiesis of Indeterminate Potential) | Mutations in the patient's own blood cells that mimic tumor signal. The #1 source of false positives in cfDNA MRD. |
| **Buffy coat** | The white-blood-cell layer of a spun blood sample. Sequencing it tells us the patient's germline AND their CHIP, so we can subtract both. |
| **Clonal / subclonal / CCF** | Clonal ("truncal") mutations are in every tumor cell; subclonal ones only in some. CCF (Cancer Cell Fraction) measures this. Truncal mutations make the best trackers. |
| **Purity / ploidy** | Tumor purity = fraction of the biopsy that is actually tumor. Ploidy = how many copies of the genome the tumor carries. Both are needed to interpret mutation fractions correctly. |
| **Allele-specific enrichment** | 2Strands' core technology: probes that selectively pull the MUTANT version of a fragment out of the blood, concentrating the signal and cutting sequencing cost. |
| **Panel** | The personalized set of mutation sites (and their probes) chosen for one patient. |
| **Genome equivalents (GE)** | How many copies of the genome are in a blood sample. ~30 ng of cfDNA ~ ~9,000 GE. This hard-caps how low a VAF you can possibly detect. |
| **PoN** (Panel of Normals) | A reference set built from healthy samples used to model background noise. There are two different ones here (Section 6.5). |
| **LoB / LoD / LoQ** | Limit of Blank / Detection / Quantification: the statistical floors that define what "negative", "detected", and "measurable" actually mean. |

---

## 6. The technologies, and what each one is for

The full pipeline is a chain of best-in-class genomics tools plus a custom
statistical engine. Standard tools do the standard jobs; the custom Python
(`bin/`) does the parts that are specific to ultra-low-frequency, enriched MRD.

| Stage | Tool | What it does | Why this one |
|---|---|---|---|
| Workflow engine | **Nextflow** (DSL2) | Orchestrates every step, parallelizes, pins container versions | Reproducible and clinical-scaling; the field standard (nf-core) |
| UMI handling + consensus | **fgbio** | Extracts UMIs, groups reads by molecule, builds duplex consensus | The reference implementation for duplex error correction |
| Alignment | **bwa-mem2** | Maps reads to the human genome (GRCh38) | Fast, deterministic when pinned; clinical-grade |
| Somatic calling (WES) | **GATK Mutect2** | Finds the tumor's mutations vs the matched normal | Mature tumor-normal somatic caller |
| Copy number | **CNVkit / GATK CNV** | Measures gains/losses of genomic regions | Feeds purity/ploidy and clonality |
| Purity / ploidy | **PURPLE / FACETS / TITAN** | Estimates tumor purity and genome copies | Needed to convert raw VAFs into cellular fractions |
| Clonality | **PyClone-vi** | Computes CCF; flags truncal vs subclonal | Picks the most reliable trackers for the panel |
| Annotation | **Ensembl VEP** | Labels each mutation (gene, effect) | Versioned, auditable annotation for a regulated setting |
| Germline filter | **gnomAD / dbSNP** | Removes common inherited variants | Population databases catch germline the matched normal might miss |
| Metrics | **Picard CollectHsMetrics** | On-target coverage QC | Standard capture-panel QC |
| Custom engine | **Python (numpy, scipy, pysam)** | Panel selection, interrogation, the MRD statistics | The MRD-specific logic that no off-the-shelf tool provides |

---

## 7. Design decisions and the reasons behind them

This is the part worth reading slowly. Each decision exists because a simpler
choice fails in a specific, concrete way.

### 7.1 Two pipelines, and interrogation instead of variant calling

A natural first instinct is to run a somatic variant caller (like Mutect2) on the
blood sample. This is wrong for the monitoring step. General somatic callers have
a practical detection floor around 1-2% VAF and built-in filters that actively
suppress anything rarer. They would throw away the 0.05% signal we are hunting
for.

Because the assay is tumor-informed, we already KNOW the patient's mutations from
Pipeline A. So Pipeline B does not "call variants", it **interrogates** known
sites: for each known mutation, count how many molecules support it and ask
whether that exceeds the site's error floor. Calling is a discovery problem;
MRD is a genotyping-at-known-sites problem. Mixing the two is why the design is
split into two pipelines (A = discovery, B = interrogation).

### 7.2 Duplex UMIs (and why paired-end is not duplex)

To suppress sequencing errors to the 1e-7 level needed for sub-0.1% detection,
the assay uses **duplex** consensus: a true mutation must appear on BOTH the
Watson and Crick strands of the original molecule. A common confusion is that
"paired-end sequencing" already does this. It does not. Paired-end means you read
both ends of one fragment; duplex means you separately tag and confirm both
complementary strands. Duplex needs specific adapter chemistry (assumed present
here). Without it, the honest detection floor is closer to 0.1-1%.

### 7.3 Allele-specific enrichment, and odds-space de-biasing

2Strands' probes selectively capture the mutant allele, which concentrates the
signal and lets the assay run "at a fraction of the sequencing burden". This is
the company's edge, but it confounds the math: the observed VAF is inflated by a
per-probe **enrichment factor** (often 10-1000x), so a raw measurement of 9% might
be a true tumor fraction of just 0.2%.

The fix models enrichment correctly: it multiplies the mutant-allele **odds**
(`odds_obs = enrich x odds_true`, where odds = vaf/(1-vaf)). The pipeline
calibrates each probe's enrichment factor from a contrived dilution series and
then de-biases reported tumor fraction in odds space. An earlier version fit the
factor in plain VAF space, which saturated at high calibration points and
under-recovered enrichment by ~2.3x; the odds model removed that bias (recovery
went from 0.43x to ~1.0x of truth, and a valid Limit of Quantification appeared).

### 7.4 Panel-integrated detection (no single site carries the call)

At 0.01% VAF, no single mutation site has enough molecules to be confident on its
own. The detection has to **integrate signal across the whole panel** of mutations
into one patient-level statistic. The engine sums the mutant-supporting molecules
across all sites and compares that total to a Monte-Carlo null distribution,
producing one p-value, one tumor-fraction estimate, and a confidence interval per
timepoint. Detection is always a composite call, never per-variant. This is what
makes sub-0.1% sensitivity possible at all, and it is bounded by the input mass:
sensitivity is limited by genome equivalents in the tube, not by the algorithm.

### 7.5 Two different backgrounds, never conflated

"What does noise look like in a sample without the tumor?" has two different
answers here, for two different pipelines:

- **Mutect2 Panel of Normals** (Pipeline A): built from normal WES samples,
 filters recurrent technical artifacts during somatic discovery.
- **Healthy-donor cfDNA empirical null** (Pipeline B): built from healthy
 blood samples run through the same assay, models the per-site error
 distribution for the MRD statistics.

They are different analytes (tissue WES vs blood cfDNA) for different purposes.
The pipeline needs both, and contrived reference samples (cell-line dilutions) are
for validation only, not as the running null.

### 7.6 The empirical, covariance-preserving null

The MRD null model is the heart of specificity (avoiding false positives). The
naive approach draws each site's error independently. But real cfDNA errors are
**correlated** across sites: a bad sequencing run, GC bias, or deamination damage
pushes many sites up together. An independent null cannot see that correlation,
so it under-estimates how often a tumor-free sample looks suspicious, and it
raises false positives.

The fix builds the null by **resampling whole healthy-donor panels**: each null
replicate draws one real donor's entire per-site error pattern as a unit, which
preserves the cross-site correlation. On the validation harness, this dropped the
false-positive rate from ~28% (independent null) back to ~7% (near the nominal 5%
target) with no loss of sensitivity. The independent null is kept only as a
back-compatible fallback that warns when used.

### 7.7 CHIP handled by biology, not by a gene list

Clonal hematopoiesis mutations in the patient's blood cells are the dominant
false-positive risk. The correct defense is to sequence the **buffy coat** (white
cells) and subtract anything present there, because CHIP is defined by presence in
the blood compartment, not by which gene it hits. (A known refinement still on the
list: today's panel-selection also drops variants by CHIP gene name, which is too
blunt because genes like DNMT3A and TET2 are both CHIP genes and real tumor
drivers. See [TODO.md](TODO.md).)

### 7.8 Probe-feasibility-aware panel selection

Not every biologically informative mutation can be turned into a good
allele-specific probe (some have bad GC content, homopolymer runs, or poor
specificity). So panel selection is the **intersection** of "informative" (truncal,
high-CCF, not CHIP, not germline) and "enrichable" (a designable, specific probe
exists). The selector pre-screens feasibility and ranks survivors by predicted
enrichment, then hands candidates to the existing 2Strands probe-design workflow.

### 7.9 Sample identity and the patient-lock (fail-closed)

A tumor-informed assay has a failure mode that has nothing to do with statistics:
running a patient's blood against **the wrong patient's panel** (a sample swap, a
mislabeled tube). The defense is two independent mechanisms, both implemented in
`bin/sample_id.py`:

- **SNP-fingerprint concordance.** Genotype a fixed common-SNP set in two BAMs and
 compare genotype dosage; a matched tumor/normal pair reads SAME, an unrelated
 pair reads DIFFERENT. The `concordance` subcommand exits non-zero on a confirmed
 DIFFERENT verdict, so a workflow step can gate on it.
- **Patient-lock provenance token.** Pipeline A's `PANEL_SELECT` stamps every panel
 with a `*.panel.lock` = `panel-<N>-<sha256>`, a deterministic hash over the
 panel's content and patient id. Pipeline B's `MRD_INTEGRATE` runs a **fail-closed**
 `sample_id.py verify-lock` *before* scoring: a panel whose lock does not match the
 cfDNA sample's patient aborts the task rather than producing a number. The matched
 token is then stamped into `mrd.json` (`mrd_integrate.py --panel-lock`) so the
 result carries proof of which panel produced it.

---

## 8. The custom engine (`bin/`)

**Eight self-tested Python scripts** hold the MRD-specific logic (the core engine
plus the glue that adapts standard-tool output to it). Each has a pure,
unit-testable core and a `selftest` subcommand that runs on synthetic data with no
sequencing required; `bash bin/run_selftests.sh` runs all eight and exits non-zero
if any fails.

Core engine:

- **`panel_select.py`** (Pipeline A) - turns annotated somatic variants into the
 personalized panel (BED + VCF). Applies SNV-only, CHIP and buffy filters,
 gnomAD/dbSNP germline exclusion, minimum-CCF, and probe feasibility; ranks by
 clonality, CCF, and predicted enrichment.
- **`interrogate.py`** (Pipeline B, stage B5) - counts unique consensus molecules
 supporting each panel mutation in a cfDNA sample, tracking duplex support and
 consensus quality. Pure counting core is BAM-free testable.
- **`mrd_integrate.py`** (Pipeline B, stage B6) - the MRD engine. Empirical
 covariance-preserving null, Monte-Carlo p-value, odds-space enrichment
 de-biasing, bootstrap confidence interval, and the Limit-of-Detection gate that
 returns "indeterminate" rather than a false "negative" when molecular depth is
 too low. Stamps the patient-lock token into `mrd.json`.

QC / identity / background:

- **`validate.py`** - the offline validation harness. Simulates blank and dilution
 cohorts (or reads real ones), fits the background on a training fold, evaluates
 a disjoint held-out fold, and reports false-positive rate, Limit of Blank,
 LoD95, LoQ, and enrichment recovery. It is deliberately built to be able to
 FAIL the engine: a clean report would mean the harness was circular.
- **`sample_id.py`** - sample-identity and patient-lock (Section 6.9): SNP-fingerprint
 genotyping + `concordance` (sample-swap detection) and `provenance`/`verify-lock`
 (the panel lock token).
- **`build_background.py`** - builds the healthy-donor cfDNA empirical null
 (per-site error distribution) that `mrd_integrate.py` consumes.

Pipeline A glue:

- **`pyclone_prep.py`** - adapts copy-number + somatic VCF into PyClone-vi input and
 converts PyClone-vi output into the per-variant CCF / clonal-probability table
 `panel_select.py` expects.
- **`normal_evidence.py`** - summarizes matched-normal / buffy-coat support per
 somatic site so `panel_select.py` can subtract germline and CHIP.

Alongside these are **real-run driver scripts** (not part of the self-tested
engine): `run_panel_design.sh`, `run_validation_chain.sh`, `fetch_url.py` (a
resumable `urllib` downloader, because `curl`/`wget` are blocked in this
environment), `fetch_wes.sh`, `queue_panel_design.sh`, and `build_panel.sh`.

See [README.md](README.md) for exact commands.

---

## 9. How sensitivity actually works (and its hard limit)

It is worth being honest about physics. Detecting 0.01% VAF means finding roughly
one mutant molecule among ten thousand. Two things make it possible:

1. **Allele-specific enrichment** concentrates the mutant molecules so they are
 over-represented in what gets sequenced.
2. **Panel integration** pools weak evidence from dozens of mutation sites into
 one call.

But there is a floor nothing can cross: the number of **genome equivalents** in
the blood tube. A standard ~30 ng draw is ~9,000 genome copies. You cannot detect
a mutant fraction so low that the expected number of mutant molecules across the
whole panel is much less than one. That is why the engine reports
**"indeterminate"** (not "negative") when molecular depth is insufficient, and why
the validation harness measures a real, input-bounded LoD rather than asserting a
marketing number.

---

## 10. Status and limitations (read before trusting any output)

This repository now contains the **custom statistical engine** (`bin/`,
self-tested) AND the **built Nextflow orchestration** for both pipelines. It is
still research-stage, not a validated clinical assay.

- **Both DSL2 pipelines are built.** `main.nf` routes
 `--workflow panel_design | mrd_monitor` to 2 workflows (`workflows/`) assembled
 from 3 subworkflows (`subworkflows/local/`) and 20 process modules
 (`modules/local/`). Containers follow the SVcaller convention: one pinned
 biocontainer per standard tool (`conf/docker.config`, tags include the conda
 build hash to avoid solver lockups and image balloon) plus one slim
 `mrd-umi/utils:1.0` image for the whole custom Python engine.
- **Engine verified:** all eight Python tools and their self-tests, plus the
 validation harness with measured before/after numbers (Sections 6.3, 6.6).
- **Sample identity / patient-lock is implemented and wired** (Section 6.9), not a
 gap any more.

### 10.1 Validation results (what has actually been run on real data)

**Pipeline A - validated end-to-end on real SEQC2 WES.** Tumor HCC1395
(SRR7890850) + matched normal HCC1395BL (SRR7890851), GRCh38, run through
align -> Mutect2 -> FACETS -> VEP -> PyClone-vi -> panel selection. Results:

- A 50-site personalized panel (committed at `assets/example_panel_HCC1395/`).
- FACETS purity **0.90** / ploidy **3.06**, matching the known near-pure aneuploid
 HCC1395 cell line - a cross-check that the run is real, not stubbed.
- Somatic-SNV **F1 = 0.79** against the SEQC2 high-confidence truth set, evaluated
 **exome-restricted** (precision 0.79 / recall 0.80). The naive genome-wide F1 of
 **0.054** is a WES-vs-whole-genome-truth artifact (WES covers ~1-2% of the
 genome, so ~38k truth SNVs are physically uncovered -> false FN); restricting the
 truth to HC regions intersect tumor depth >=10x (from `mosdepth`) is what makes the
 comparison fair. Full method and numbers live in
 [VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md).

**Pipeline B - plumbing proven, LoD still pending.** A downsampled
(`fastq-dump -X`) SEQC2 ILM2 titration chain proves the cfDNA
download -> align -> interrogate path on real reads, but per-site coverage is too
sparse for a real limit-of-detection. Full-depth LoD/LoQ remains pending (it is
egress-bound on this box).

### 10.2 A lesson: `-stub-run` validates wiring, not behavior

The Nextflow DAG was green under `-stub-run` (stub blocks just `touch` their
outputs, so they prove channel/join wiring independent of tool behavior or
reference contents) long before the first real run. A series of bugs that stubs
**cannot** catch only surfaced on the real run, and were fixed:

1. Reference companions not staged with the FASTA (bwa-mem2 index, the SNP `.tbi`,
 the `.fai`).
2. bwa-mem2 `-C` corrupting WES SAMs - `-C` is only valid for the interleaved
 UMI/RX flow, not paired WES.
3. `pyclone_prep` defaulting tumor to `samples[0]`, which Mutect2 had ordered as
 the NORMAL (-> CCF 0 -> empty panel); fixed by passing `--tumor-sample`.
4. A non-executable `bin/` script (exit 126).
5. Mutect2 PoN / germline-resource `.tbi` indices not staged.

The takeaway: stub-runs are necessary but not sufficient; only a real run exercises
reference staging, tool flags, and sample ordering.

### 10.3 Project plumbing

- **CI** (`.github/workflows/ci.yml`) runs the self-tests and both-pipeline
 `-stub-run` DAG checks - data-independent, no containers required.
- A **PostToolUse hook** (`bin/hooks/selftest_on_edit.sh`) runs the matching
 self-test on each `bin/*.py` edit, so a regression surfaces in-loop.
- **`publishDir`** copies the panel to `results/panel_design` and `mrd.json` to
 `results/mrd`.

### 10.4 Still required before clinical use

- **Important caveat on the self-tests:** they prove the code correctly implements
 its model, not that the model matches real biology. Real validity comes only
 from running `validate.py run-real` on wet-lab dilution series and healthy-donor
 cohorts, and from a full-depth cfDNA LoD.
- **Not yet present:** run-level controls, the CHIP gene-list refinement
 (Section 7.7), fragment-end trimming, a complete provenance/audit block (the
 patient-lock is stamped, but reference hash / container digests / code version /
 RNG seed are not yet), and SOPs. These are tracked in [TODO.md](TODO.md).

Nothing here is for clinical decision-making in its current form.

## References (selected methodology)

Datasets used for validation are cited in
[VALIDATION_DATASETS.md](VALIDATION_DATASETS.md). The methods this pipeline builds
on:

**Tumor-informed ctDNA / MRD**
- Abbosh C, et al. Phylogenetic ctDNA analysis depicts early-stage lung cancer evolution. *Nature*. 2017;545:446-451.
- Wan JCM, et al. Liquid biopsies come of age: towards implementation of circulating tumour DNA. *Nat Rev Cancer*. 2017;17:223-238.

**UMI / duplex consensus (error correction)**
- Schmitt MW, et al. Detection of ultra-rare mutations by next-generation sequencing. *Proc Natl Acad Sci USA*. 2012;109:14508-14513.
- Kennedy SR, et al. Detecting ultralow-frequency mutations by Duplex Sequencing. *Nat Protoc*. 2014;9:2586-2606.
- fgbio. Fulcrum Genomics. https://github.com/fulcrumgenomics/fgbio

**Background / position-specific error model (the empirical null)**
- Newman AM, et al. An ultrasensitive method for quantitating circulating tumor DNA with broad patient coverage. *Nat Med*. 2014;20:548-554.
- Newman AM, et al. Integrated digital error suppression for improved detection of circulating tumor DNA. *Nat Biotechnol*. 2016;34:547-555.

**Somatic discovery, copy number, clonality, annotation (Pipeline A)**
- Benjamin D, et al. Calling Somatic SNVs and Indels with Mutect2. *bioRxiv*. 2019. doi:10.1101/861054.
- Shen R, Seshan VE. FACETS: allele-specific copy number and clonal heterogeneity analysis tool. *Nucleic Acids Res*. 2016;44:e131.
- Roth A, et al. PyClone: statistical inference of clonal population structure in cancer. *Nat Methods*. 2014;11:396-398.
- Gillis S, Roth A. PyClone-VI: scalable inference of clonal population structure using variational inference. *BMC Bioinformatics*. 2020;21:571.
- McLaren W, et al. The Ensembl Variant Effect Predictor. *Genome Biol*. 2016;17:122.
- Vasimuddin M, et al. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *IEEE IPDPS*. 2019:314-324.

**Clonal hematopoiesis (CHIP filtering)**
- Jaiswal S, et al. Age-related clonal hematopoiesis associated with adverse outcomes. *N Engl J Med*. 2014;371:2488-2498.

**Allele-specific / minor-allele enrichment (the enrichment layer's method class)**
- Song C, et al. Elimination of unaltered DNA in mixed clinical samples via nuclease-assisted minor-allele enrichment (NaME-PrO). *Nucleic Acids Res*. 2016;44:e146.
- Li J, et al. Replacing PCR with COLD-PCR enriches variant DNA sequences and redefines the sensitivity of genetic testing. *Nat Med*. 2008;14:579-583.

**Analytical validation (limits of blank/detection/quantitation)**
- CLSI. Evaluation of Detection Capability for Clinical Laboratory Measurement Procedures. Guideline EP17-A2. Clinical and Laboratory Standards Institute; 2012.
