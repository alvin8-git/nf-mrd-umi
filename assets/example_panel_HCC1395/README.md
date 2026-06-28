# Example panel - HCC1395 (real Pipeline A output)

Validation evidence: the personalized MRD panel produced by **Pipeline A
(`panel_design`)** run end-to-end on real SEQC2 WES - tumor **HCC1395**
(SRR7890850) + matched normal **HCC1395BL** (SRR7890851), GRCh38, with a Mutect2
**1000G panel-of-normals + gnomAD** germline resource.

- `HCC1395.panel.bed` - 50 selected sites (`chrom  start  end  GENE:REF>ALT`)
- `HCC1395.panel.vcf` - the same 50 variants, VEP-annotated
- `HCC1395.panel.lock` - patient-lock provenance token (`panel-50-a6aef193fab5e1db`);
  Pipeline B verifies this before interrogating a cfDNA sample (sample-swap protection)

Pipeline cross-checks that this run is real, not stubbed:
- FACETS purity **0.90** / ploidy **3.06** - matches the known near-pure, aneuploid
  HCC1395 cell line.
- **50 / 50** selected panel sites are confirmed SEQC2 high-confidence somatic SNVs
  (up from 33/50 before the PoN+gnomAD layer - the germline/artifact contaminants
  are gone, so every tracked site is a bona-fide tumor mutation).
- Somatic-calling F1 = **0.806** in the exome-callable HC region (precision 0.822 /
  recall 0.790); see [VALIDATION_COVERAGE.md](../../VALIDATION_COVERAGE.md).

Regenerate with `bin/run_panel_design.sh` (set `PON=` + `GERMLINE=`; publishes to
`results/panel_design/`).
