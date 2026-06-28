#!/usr/bin/env bash
# fetch_wes.sh - download SEQC2 HCC1395/HCC1395BL WES runs (full depth) for Pipeline A.
# prefetch -> verify LOCAL .sra (else fasterq-dump streams remotely) -> paired fastq.
# Usage: PATH=/data/alvin/envs/mrd/bin:$PATH bash bin/fetch_wes.sh OUTDIR SRR [SRR...]
set -euo pipefail
out=${1:?usage: fetch_wes.sh OUTDIR SRR [SRR...]}; shift
threads=${THREADS:-8}
mkdir -p "$out"
for srr in "$@"; do
  if [[ -s "$out/${srr}_1.fastq" || -s "$out/${srr}_1.fastq.gz" ]]; then
    echo "[$srr] skip (fastq exists)"; continue
  fi
  echo "[$srr] prefetch -> local .sra"
  prefetch -O "$out" --max-size 100G "$srr" >/dev/null
  sra=$(find "$out" -name "${srr}.sra" -o -name "${srr}.sralite" 2>/dev/null | head -1)
  [[ -s $sra ]] || { echo "[$srr] ERROR: no local .sra (would stream)" >&2; exit 1; }
  echo "[$srr] extract $(du -h "$sra"|cut -f1) -> fastq"
  mkdir -p "$out/fqtmp.$srr"
  fasterq-dump --split-files --threads "$threads" -O "$out" -t "$out/fqtmp.$srr" "$sra" >/dev/null
  rm -rf "$out/fqtmp.$srr"
  [[ -s "$out/${srr}_1.fastq" ]] || { echo "[$srr] ERROR: no fastq produced" >&2; exit 1; }
  rm -f "$sra"   # reclaim space; fastq is what Pipeline A consumes
  echo "[$srr] done -> $out/${srr}_{1,2}.fastq"
done
echo "all WES runs fetched -> $out"
