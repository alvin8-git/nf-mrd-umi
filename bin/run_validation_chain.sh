#!/usr/bin/env bash
# run_validation_chain.sh — headless: full downsampled SEQC2 batch -> background -> validate.
# Waits for any in-flight single-run smoke (so workdirs don't collide), then drives the
# whole manifest (the driver resumes, skipping runs whose output already exists),
# builds background.tsv/pon.tsv from the `blank` rows, and runs validate.py run-real.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

export PATH=/data/alvin/envs/mrd/bin:$PATH
MANIFEST=data/manifests/seqc2_ilm2_titration.manifest.tsv
REF=/data/alvin/ref/GRCh38/hg38.fa
PANEL=panel.vcf.gz
WORK=/data/alvin/tmp/smoke_work
OUT=results/interrogate
N=3000000

# 1) wait out any running smoke driver (different process; poll, not `wait`)
for pid in $(pgrep -f 'fetch_and_interrogate.sh' || true); do
  [[ $pid == $$ ]] && continue
  echo "[chain] waiting for in-flight driver PID $pid"
  while kill -0 "$pid" 2>/dev/null; do sleep 30; done
done

# 2) full downsampled batch (resumes; skips runs with existing output)
echo "[chain] downsampled batch over $(wc -l <"$MANIFEST") manifest rows"
bash bin/fetch_and_interrogate.sh -m "$MANIFEST" -r "$REF" -p "$PANEL" \
  -t 16 -N "$N" -w "$WORK" -o "$OUT"

# 3) background + PoN from the blank rows
echo "[chain] build background from blank rows"
python3 bin/build_background.py run --manifest "$MANIFEST" --label blank \
  --out-background background.tsv --out-pon pon.tsv

# 4) validate on real reads
echo "[chain] validate.py run-real"
python3 bin/validate.py run-real --manifest "$MANIFEST" \
  --background background.tsv --pon pon.tsv \
  --out results/seqc2_validation.json

echo "[chain] DONE -> results/seqc2_validation.json"
