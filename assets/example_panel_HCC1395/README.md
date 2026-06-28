# Example panel — HCC1395 (real Pipeline A output)

Validation evidence: the personalized MRD panel produced by **Pipeline A
(`panel_design`)** run end-to-end on real SEQC2 WES — tumor **HCC1395**
(SRR7890850) + matched normal **HCC1395BL** (SRR7890851), GRCh38.

- `HCC1395.panel.bed` — 50 selected sites (`chrom  start  end  GENE:REF>ALT`)
- `HCC1395.panel.vcf` — the same 50 variants, VEP-annotated
- `HCC1395.panel.lock` — patient-lock provenance token (`panel-50-...`); Pipeline B
  verifies this before interrogating a cfDNA sample (sample-swap protection)

Pipeline cross-checks that this run is real, not stubbed:
- FACETS purity **0.90** / ploidy **3.06** — matches the known near-pure, aneuploid
  HCC1395 cell line.
- 50 clonal SNVs survived CHIP / buffy / gnomAD / CCF filtering (min_ccf 0.10).

Regenerate with `bin/run_panel_design.sh` (now publishes to `results/panel_design/`).
See [[../../VALIDATION_COVERAGE.md]] for the full component->dataset map.
