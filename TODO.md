# nf-mrd-umi: Roadmap & TODO

Status snapshot and prioritized work. Priorities come from an independent 3-lens
design review (statistics, assay biology, clinical/regulatory). Background for
every item is in [Documentation.md](Documentation.md).

Legend: `[x]` done | `[~]` partial | `[ ]` not started

---

## Done

- [x] Two-pipeline architecture designed (A: panel design, B: MRD monitoring)
- [x] `panel_select.py` - personalized panel selection (SNV-only, CHIP + buffy +
      gnomAD/dbSNP germline filters, probe-feasibility-aware, clonality ranking)
- [x] `interrogate.py` - duplex-aware unique-molecule counting at panel sites
- [x] `mrd_integrate.py` - MRD engine: empirical null, Monte-Carlo p-value,
      odds-space enrichment de-biasing, bootstrap CI, LoD/indeterminate gate
- [x] `validate.py` - offline validation harness (simulate-run, run-real, selftest)
- [x] **Empirical covariance-preserving null** - fixed the anticonservative
      independent null (validated false-positive rate 0.28 -> 0.07 at nominal 0.05)
- [x] **Odds-space enrichment calibration** - fixed the saturating VAF-space fit
      (enrichment recovery 0.43x -> 1.0x; valid Limit of Quantification appears)
- [x] Self-tests for all eight engine scripts (`bin/run_selftests.sh`)
- [x] **Both Nextflow DSL2 pipelines built** - `panel_design` (A) + `mrd_monitor`
      (B): 20 modules, 3 subworkflows, 2 workflows, `-stub-run` DAG-validated
- [x] **Pipeline A validated end-to-end on real SEQC2 HCC1395/HCC1395BL WES** -
      50-site panel; FACETS purity 0.90 / ploidy 3.06 (matches truth);
      somatic-SNV **F1 = 0.79** exome-restricted (precision 0.79 / recall 0.80).
      Method + numbers in [VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md); the
      validated panel is committed at `assets/example_panel_HCC1395/`.
- [x] **Six real-run bugs fixed that stub-runs cannot catch** - reference-companion
      staging (bwa index / `.tbi` / `.fai`), bwa-mem2 `-C` corrupting WES SAMs,
      pyclone_prep picking the normal sample, exec-bit, Mutect2 PoN/germline `.tbi`.
- [x] CI (`.github/workflows/ci.yml`: selftests + both-pipeline stub-runs),
      PostToolUse selftest hook, `nf-pipeline-reviewer` subagent, `publishDir`.

---

## P0 - Clinical defensibility (do before any real patient sample)

- [x] **Sample-identity concordance.** `bin/sample_id.py` - SNP-fingerprint
      genotyping + concordance (`fingerprint`/`concordance`, exits non-zero on a
      confirmed DIFFERENT verdict so a step can gate). Validates against matched
      vs unrelated samples.
- [x] **Patient-locked panel artifact + verify.** `PANEL_SELECT` emits a
      content+patient hash (`*.panel.lock`, `panel-<N>-<sha>`); `MRD_INTEGRATE`
      runs a fail-closed `sample_id.py verify-lock` before scoring (refuses a
      mismatched panel / wrong patient).
- [x] **Provenance stamp in `mrd.json`.** The panel-lock token is written into the
      result (`mrd_integrate.py --panel-lock`). _Still to add: reference hash,
      container digests, code version, RNG seed into the same block._
- [~] **Validation harness on real data.** Pipeline A done (somatic F1 above).
      Pipeline B: downsampled SEQC2 ILM2 titration chain (`run_validation_chain.sh`)
      proves the cfDNA plumbing on real reads; **full-depth LoD/LoQ still pending**
      (egress-bound; downsampled coverage is too sparse for a real limit). Datasets
      in [VALIDATION_DATASETS.md](VALIDATION_DATASETS.md).
- [ ] **Run-level controls.** Positive control, no-template control, and a
      contrived-VAF reference per run, with pass/fail QC gating that fails the
      whole run (not just one sample).

---

## P1 - Specificity & correctness hardening

- [ ] **CHIP gene-blocklist -> flag, not drop.** Genes like DNMT3A/TET2 are both
      CHIP genes and real tumor drivers; the current hard gene-drop discards
      genuine tumor variants. Rely on buffy subtraction; downgrade the gene list
      to an annotation flag.
- [ ] **Lower the buffy-coat normal VAF floor.** Today's 2% floor lets 0.5-2% CHIP
      clones leak into the 0.01-0.1% detection band. Drop to ~0.1-0.5% with
      error-corrected, matched-depth buffy.
- [ ] **End-repair artifact trimming.** Trim ~5-10 bp from cfDNA fragment ends in
      `interrogate.py`; end-repair fill-in can create TRUE duplex-supported false
      positives that duplex consensus cannot remove.
- [ ] **Use duplex support in the call.** `interrogate.py` records `duplex_support`
      but `mrd_integrate` pools simplex + duplex; model and weight them separately.
      Verify fgbio `aD`/`bD` tag semantics (per-base arrays, not scalars).
- [ ] **Sample-specific LoD + calibrated alpha.** Replace the fixed molecule floor
      with a per-sample LoD computed from that sample's depths/backgrounds; set
      alpha from the empirical Limit of Blank (a flat 0.05 per timepoint compounds
      to ~26% over six serial draws).
- [ ] **Robust detection statistic.** Consider a per-site likelihood-ratio or
      winsorized contribution so one noisy high-depth site cannot dominate; require
      >=2 informative sites.
- [ ] **`_gene_from_csq` parse by field name**, not position (VEP CSQ order varies;
      a positional parse silently mis-assigns genes and breaks the CHIP filter).

---

## P2 - Sensitivity & scope

- [ ] **Minimum-trackable-N gate + tumor-type eligibility.** Refuse to ship an
      underpowered panel; flag low-mutation-burden, CNA-driven, and fusion-driven
      tumors (e.g. RCC, thyroid, sarcoma) as assay-ineligible rather than silently
      shipping too few trackers.
- [ ] **Fragmentomics capture.** Record fragment length (TLEN) per molecule in
      `interrogate.py` now (cheap, future-proofing); add a short-fragment weighting
      model later (tumor cfDNA is shorter, ~134-144 bp vs ~167 bp).
- [ ] **Input-mass / genome-equivalents gate.** Gate on input GE and conversion
      efficiency; raise the molecular-depth floor to a level consistent with the
      sub-0.1% claim.
- [ ] **Indels in the panel (v2).** Add high-confidence somatic indels with their
      own error model for patients with few trackable SNVs.

---

## Build-out (the orchestration layer)

- [x] **Nextflow DSL2 implementation.** `main.nf` routes `--workflow
      panel_design | mrd_monitor`; `workflows/`, `subworkflows/local/`,
      `modules/local/` (nf-core style); one pinned biocontainer per tool in
      `conf/docker.config` + one slim `mrd-umi/utils:1.0` engine image.
- [x] **`--intervals` (exome BED) wired** into Mutect2 (optional). Passing it cuts
      the whole-genome Mutect2 scan (~90 min) to minutes and sharpens WES calling.
- [ ] **`probe_design.nf` subworkflow.** Interface to (adapt) the existing 2Strands
      probe-design workflow; feed designable/specific/enrichment scores back into
      `panel_select`. (`panel_select.py` already consumes probe scores; the
      upstream design step is external.)
- [~] **Reference + determinism lock.** Containers are digest-pinnable; still to do:
      pin GRCh38 by hash and pin bwa-mem2 threading (thread-count non-determinism).
- [ ] **`nf-test` fixtures** wrapping the modules, including a tiny BAM fixture so
      `interrogate.py`'s pysam path gets exercised (today only its pure core is,
      via the CI stub-run + selftests).
- [ ] **Mutect2 best-practice resources by default.** A re-run with a 1000G PoN +
      gnomAD germline resource is in progress to lift somatic F1 from 0.79 toward
      ~0.9 (the gap is germline/artifact FPs a PoN+gnomAD remove).

---

## External dependencies (not code we write)

- [ ] Healthy-donor cfDNA cohort for the empirical null (sized; disjoint from the
      specificity-validation fold).
- [ ] Contrived VAF dilution series for enrichment calibration and LoD/LoQ.
- [ ] Confirmation that the enrichment chemistry is duplex-capable.
- [ ] SOPs and validation documentation for the regulated (CAP/CLIA) setting.

---

## Open scientific question

- [ ] **Validate the odds enrichment model against real chemistry.** The pipeline
      assumes allele-specific enrichment acts as a mutant-allele odds multiplier.
      The simulation cannot confirm this; only `validate.py run-real` on real
      contrived dilutions can. If the real relationship differs, the de-biasing
      model must be refit.
