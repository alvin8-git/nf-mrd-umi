#!/usr/bin/env bash
# queue_panel_design.sh - headless: block until the VEP cache is ready, then
# pause the ILM2 chain (resumable), run Pipeline A on the HCC1395 WES pair, and
# resume the ILM2 chain afterward. Per the agreed plan: Pipeline A is a one-shot
# run that unlocks 3 validations, so it briefly preempts the ILM2 LoD sweep.
#
# This script's name matches NONE of the ILM2 pgrep patterns, so its pkill calls
# below cannot match itself.
set -uo pipefail
export PATH=/data/alvin/envs/mrd/bin:/home/alvin/bin:$PATH
cd "$(dirname "${BASH_SOURCE[0]}")/.."
log() { echo "[queue $(date +%H:%M:%S)] $*"; }

FETCH_LOG=/data/alvin/tmp/resource_fetch.log
log "waiting for VEP_CACHE_READY in $FETCH_LOG"
until grep -q VEP_CACHE_READY "$FETCH_LOG" 2>/dev/null; do
  grep -q INCOMPLETE "$FETCH_LOG" 2>/dev/null && { log "ABORT: resource fetch reported INCOMPLETE"; exit 1; }
  sleep 60
done
log "VEP cache ready"

# pause the ILM2 chain (it resumes later, skipping completed tsvs)
pkill -9 -f 'run_validation_chain.sh' 2>/dev/null || true
pkill -9 -f 'fetch_and_interrogate.sh' 2>/dev/null || true
pkill -9 -f 'fastq-dump -X 3000000' 2>/dev/null || true
log "paused ILM2 chain ($(ls results/interrogate/*.tsv 2>/dev/null | wc -l)/48 tsvs banked)"

# launch Pipeline A (real run; first use pulls gatk4/picard/vep/facets/pyclone-vi)
export SNP_VCF=/data/alvin/ref/GRCh38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
export VEP_CACHE=/data/alvin/ref/vep_cache
# best-practice Mutect2 resources (improve precision); used if present
[[ -f /data/alvin/ref/GRCh38/af-only-gnomad.hg38.vcf.gz ]] && export GERMLINE=/data/alvin/ref/GRCh38/af-only-gnomad.hg38.vcf.gz
[[ -f /data/alvin/ref/GRCh38/1000g_pon.hg38.vcf.gz ]] && export PON=/data/alvin/ref/GRCh38/1000g_pon.hg38.vcf.gz
log "launching Pipeline A (panel_design) -> results/panel_design"
bash bin/run_panel_design.sh; rc=$?
log "Pipeline A exited rc=$rc"

# resume the ILM2 chain regardless of Pipeline A outcome
log "resuming ILM2 chain"
nohup bash bin/run_validation_chain.sh > /data/alvin/tmp/validation_chain.log 2>&1 &
log "ILM2 chain resumed (pid $!); queue done"
