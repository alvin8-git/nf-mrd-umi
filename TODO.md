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
- [x] Self-tests for all four scripts

---

## P0 - Clinical defensibility (do before any real patient sample)

- [ ] **Sample-identity concordance.** SNP-fingerprint module comparing the cfDNA
      consensus BAM against the tumor/buffy WES; hard-fail the run on mismatch.
      Without it, a sample swap produces a confident, wrong, signed-out result.
      (Biggest accreditation risk.)
- [ ] **Patient-locked, immutable panel artifact.** Hash + patient-stamp the
      BED/VCF panel; `mrd_integrate` must refuse if the panel's patient_id differs
      from `--patient-id`, and record the panel SHA256 in the result.
- [ ] **Provenance in `mrd.json`.** Write background hash, panel hash, reference
      hash, container digests, code version, and the RNG seed into the result so
      it is reproducible and defensible.
- [ ] **Validation harness on real data.** Run `validate.py run-real` against a
      held-out healthy-donor blank cohort (-> empirical LoB/specificity) and a
      contrived VAF dilution series (-> empirical LoD/LoQ). The simulator only
      proves the implementation; real data proves the model.
      **Candidate public datasets are catalogued in
      [VALIDATION_DATASETS.md](VALIDATION_DATASETS.md)** — start with SEQC2
      liquid biopsy (BioProject PRJNA677999).
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

- [ ] **Nextflow DSL2 implementation.** Turn the Part-2 module design into real
      `main.nf`, `workflows/`, `subworkflows/`, and `modules/` (nf-core style),
      with two entry points (`panel_design`, `mrd_monitor`) and digest-pinned
      containers.
- [ ] **`probe_design.nf` subworkflow.** Interface to (adapt) the existing 2Strands
      probe-design workflow; feed designable/specific/enrichment scores back into
      `panel_select`.
- [ ] **Reference + determinism lock.** Pin GRCh38 no-alt + decoy by hash; pin
      bwa-mem2 version and threading (thread-count non-determinism is a known trap).
- [ ] **`nf-test` fixtures** wrapping the modules, including a tiny BAM fixture so
      `interrogate.py`'s pysam path gets exercised (today only its pure core is).

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
