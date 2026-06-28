#!/usr/bin/env bash
# run_panel_design.sh - drafted launch for Pipeline A (panel_design) on the
# SEQC2 HCC1395 (tumor) / HCC1395BL (normal) WES pair.
#
# STAGED, not turnkey: three resource files must be provided first (see PREREQS).
# Fill the *_VCF / *_DIR paths below (or export them), then run this script.
# It fails loud if a required input is missing, so nothing silently no-ops.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

# --- present on this box (verified) ---
FASTA=/data/alvin/ref/GRCh38/hg38.fa
FAI=${FASTA}.fai
DICT=/data/alvin/ref/GRCh38/hg38.dict
SAMPLESHEET=assets/samplesheet_wes_hcc1395.csv
CHIP=assets/chip_blocklist.txt

# --- REQUIRED, must be fetched (see PREREQS at bottom). Override via env. ---
SNP_VCF=${SNP_VCF:-/data/alvin/ref/GRCh38/dbsnp_common.vcf.gz}     # FACETS pileup
VEP_CACHE=${VEP_CACHE:-/data/alvin/ref/vep_cache}                  # ensembl-vep offline cache (GRCh38)

# --- OPTIONAL (improve somatic calling; empty = skipped). Override via env. ---
GERMLINE=${GERMLINE:-}        # af-only-gnomad.hg38.vcf.gz  (Mutect2 germline resource)
PON=${PON:-}                  # 1000g_pon.hg38.vcf.gz       (Mutect2 panel-of-normals, WES)
INTERVALS=${INTERVALS:-}      # WES target BED (LL center exome kit; speeds/sharpens Mutect2)

OUTDIR=${OUTDIR:-results/panel_design}

# --- fail loud on missing required inputs ---
for f in "$FASTA" "$FAI" "$DICT" "$SAMPLESHEET" "$CHIP" "$SNP_VCF" "$VEP_CACHE"; do
  [[ -e $f ]] || { echo "MISSING required input: $f  (see PREREQS in $(basename "$0"))" >&2; exit 2; }
done

opt=()
[[ -n $GERMLINE  ]] && opt+=(--germline_resource "$GERMLINE")
[[ -n $PON       ]] && opt+=(--mutect2_pon "$PON")
[[ -n $INTERVALS ]] && opt+=(--intervals "$INTERVALS")

set -x
nextflow run main.nf -profile docker -resume \
  --workflow panel_design \
  --input "$SAMPLESHEET" \
  --fasta "$FASTA" --fasta_fai "$FAI" --fasta_dict "$DICT" \
  --snp_vcf "$SNP_VCF" \
  --chip_blocklist "$CHIP" \
  --vep_cache "$VEP_CACHE" \
  --vep_assembly GRCh38 \
  --panel_size 50 --min_ccf 0.10 \
  --outdir "$OUTDIR" \
  "${opt[@]}"

# ============================ PREREQS (fetch first) ============================
# Egress is throttled (~1-2 MB/s). VEP cache is the big one (~15 GB).
#
# 1) SNP_VCF (FACETS) - dbSNP common SNPs, GRCh38, bgzipped+tabixed. e.g.
#      GATK bundle Homo_sapiens_assembly38.dbsnp138.vcf, or dbSNP 00-common_all.vcf.gz
#      (must match hg38 'chr'-prefixed contigs, as our reference is UCSC-style).
# 2) VEP_CACHE - ensembl-vep offline cache for GRCh38:
#      vep_install -a cf -s homo_sapiens -y GRCh38 -c /data/alvin/ref/vep_cache
# 3) GERMLINE (optional) - GATK af-only-gnomad.hg38.vcf.gz
#    PON      (optional) - GATK 1000g_pon.hg38.vcf.gz
#    INTERVALS(optional) - the LL-center exome target BED (or any GRCh38 exome kit BED)
#
# NOTE: first run pulls the gatk4 / picard / ensembl-vep / cnv_facets / pyclone-vi
# biocontainers (conf/docker.config pins). bwa-mem2, samtools, utils:1.0 already present.
# Smoke the DAG without data first:  nextflow run main.nf -profile docker,test -stub-run --workflow panel_design
