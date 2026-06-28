#!/usr/bin/env bash
# fetch_and_interrogate.sh
#
# Fill a run-real manifest's site_counts_path files. For each SRA run in the
# manifest: download -> align -> [UMI duplex consensus] -> interrogate.py.
# Resumable (skips runs whose output already exists) and parallelizable (-j).
#
# This bridges a public dataset (e.g. SEQC2 PRJNA677999) to validate.py run-real.
# See data/manifests/README.md for the cohort and caveats.
#
# Prerequisites you provide:
#   - sra-tools (prefetch, fasterq-dump), bwa-mem2, samtools, python3 (numpy/scipy/pysam)
#   - for --umi: fgbio (+ Java)
#   - a bwa-mem2-indexed GRCh38 reference (-r)
#   - the tumor-informed panel VCF (-p): SEQC2 high-confidence HCC1395 calls
#     INTERSECT the assay's target BED (build this separately; it is the
#     "personalized panel" for this cell-line cohort)

set -euo pipefail

BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage: $(basename "$0") -m MANIFEST -r REF.fa -p PANEL.vcf [options]

Required:
  -m FILE   run-real manifest TSV (column 1 = SRA run accession)
  -r FILE   bwa-mem2-indexed reference FASTA (GRCh38 no-alt + decoy)
  -p FILE   panel VCF (sites to interrogate)

Options:
  -o DIR    output dir for site-count TSVs   (default: results/interrogate)
  -w DIR    scratch/work dir                 (default: work)
  -j N      parallel runs                    (default: 1)
  -t N      threads per run (bwa/samtools)    (default: 8)
  -u        UMI duplex consensus path (fgbio); requires -s
  -s STR    fgbio read structure for -u, e.g. "8M+T 8M+T" (assay-specific)
  -k        keep per-run intermediates (default: delete to save disk)
  -h        this help

Output paths are derived from the manifest's site_counts_path column shape:
  <outdir>/<RUN>.tsv   (one per run; matches results/interrogate/<RUN>.tsv)
EOF
}

OUTDIR="results/interrogate"; WORKDIR="work"; JOBS=1; THREADS=8
UMI=0; READ_STRUCT=""; KEEP=0; MANIFEST=""; REF=""; PANEL=""; MAXREADS=""
while getopts "m:r:p:o:w:j:t:us:kN:h" opt; do
  case "$opt" in
    m) MANIFEST=$OPTARG ;; r) REF=$OPTARG ;; p) PANEL=$OPTARG ;;
    o) OUTDIR=$OPTARG ;; w) WORKDIR=$OPTARG ;; j) JOBS=$OPTARG ;;
    t) THREADS=$OPTARG ;; u) UMI=1 ;; s) READ_STRUCT=$OPTARG ;;
    k) KEEP=1 ;; N) MAXREADS=$OPTARG ;; h) usage; exit 0 ;; *) usage; exit 2 ;;
  esac
done

# --- validate inputs (fail loud at the trust boundary) ---
[[ -n $MANIFEST && -f $MANIFEST ]] || { echo "ERROR: -m manifest missing" >&2; exit 2; }
[[ -n $REF && -f $REF ]]           || { echo "ERROR: -r reference missing" >&2; exit 2; }
[[ -f ${REF}.bwt.2bit.64 || -f ${REF}.0123 ]] || \
  echo "WARN: ${REF} does not look bwa-mem2-indexed (run: bwa-mem2 index $REF)" >&2
[[ -n $PANEL && -f $PANEL ]]       || { echo "ERROR: -p panel VCF missing" >&2; exit 2; }
if [[ $UMI -eq 1 && -z $READ_STRUCT ]]; then
  echo "ERROR: -u (UMI consensus) requires -s READ_STRUCTURE (e.g. '8M+T 8M+T')." >&2
  echo "       The read structure is assay-specific; do not guess it." >&2
  exit 2
fi

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing tool: $1" >&2; exit 2; }; }
need prefetch; need fasterq-dump; need bwa-mem2; need samtools; need python3
[[ $UMI -eq 1 ]] && need fgbio

mkdir -p "$OUTDIR" "$WORKDIR"
FAILLOG="$WORKDIR/failures.log"; : > "$FAILLOG"

# --- process one SRA run end to end ---
process_one() {
  local srr=$1
  local out="$OUTDIR/${srr}.tsv"
  if [[ -s $out ]]; then echo "[$srr] skip (output exists)"; return 0; fi
  local wd="$WORKDIR/$srr"; mkdir -p "$wd"
  local rg="@RG\tID:${srr}\tSM:${srr}\tLB:${srr}\tPL:ILLUMINA"
  if [[ -n $MAXREADS ]]; then
    # Downsampled mode: fastq-dump -X streams and STOPS after MAXREADS spots, so
    # the download is bounded to ~the first portion of the run (minutes, not the
    # full 40GB). Validates the pipeline on real reads at reduced depth; the
    # full-depth LoD needs the complete run on a high-egress box.
    echo "[$srr] downsampled download: first ${MAXREADS} spots (fastq-dump -X)"
    fastq-dump -X "$MAXREADS" --split-files -O "$wd" "$srr" >/dev/null
  else
    echo "[$srr] download (prefetch -> local .sra)"
    prefetch -O "$wd" --max-size 100G "$srr" >/dev/null
    # Verify a LOCAL .sra exists, else fasterq-dump silently streams remotely.
    local sra
    sra=$(find "$wd" -name "${srr}.sra" -o -name "${srr}.sralite" 2>/dev/null | head -1)
    if [[ ! -s $sra ]]; then
      echo "[$srr] ERROR: prefetch produced no local .sra (would stream remotely)" >&2
      return 1
    fi
    echo "[$srr] extract local .sra ($(du -h "$sra" | cut -f1)) -> fastq"
    mkdir -p "$wd/fqtmp"
    fasterq-dump --split-files --threads "$THREADS" -O "$wd" -t "$wd/fqtmp" "$sra" >/dev/null
    rm -rf "$wd/fqtmp"
  fi

  # locate fastqs (paired or single)
  local fq1="$wd/${srr}_1.fastq" fq2="$wd/${srr}_2.fastq" fqs
  if [[ -f $fq1 && -f $fq2 ]]; then fqs="$fq1 $fq2"
  elif [[ -f "$wd/${srr}.fastq" ]]; then fqs="$wd/${srr}.fastq"
  else echo "[$srr] ERROR: no fastq produced" >&2; return 1; fi

  local bam="$wd/${srr}.bam"
  if [[ $UMI -eq 1 ]]; then
    echo "[$srr] UMI duplex consensus (read-structure: $READ_STRUCT)"
    # 1) uBAM with UMI extracted into the RX tag
    fgbio FastqToBam --input $fqs --read-structures $READ_STRUCT \
      --sample "$srr" --library "$srr" --output "$wd/unmapped.bam"
    # 2) align (carry RX through: -T RX on fastq, -C on bwa copies it to tags)
    samtools fastq -T RX "$wd/unmapped.bam" \
      | bwa-mem2 mem -t "$THREADS" -p -C -R "$rg" "$REF" - \
      | samtools sort -@ "$THREADS" -o "$wd/mapped.bam" -
    # 3) group by UMI + position, build duplex consensus, filter
    fgbio GroupReadsByUmi --input "$wd/mapped.bam" --strategy paired \
      --output "$wd/grouped.bam"
    fgbio CallDuplexConsensusReads --input "$wd/grouped.bam" \
      --output "$wd/cons.unmapped.bam" --min-reads 1
    # ponytail: --min-reads 2 1 1 and FilterConsensusReads thresholds are
    # assay-specific QC knobs; tune for the real assay before clinical use.
    fgbio FilterConsensusReads --input "$wd/cons.unmapped.bam" --ref "$REF" \
      --output "$wd/cons.filt.bam" --min-reads 1 --min-base-quality 30
    # 4) realign the consensus reads -> final molecule-level BAM
    samtools fastq "$wd/cons.filt.bam" \
      | bwa-mem2 mem -t "$THREADS" -p -R "$rg" "$REF" - \
      | samtools sort -@ "$THREADS" -o "$bam" -
  else
    echo "[$srr] align + markdup (no-UMI path)"
    # bwa-mem2 emits mates adjacent, so fixmate works without a name sort.
    bwa-mem2 mem -t "$THREADS" -R "$rg" "$REF" $fqs \
      | samtools fixmate -m -u - - \
      | samtools sort -u -@ "$THREADS" - \
      | samtools markdup -@ "$THREADS" - "$bam"
    # CAVEAT: interrogate.py counts one read = one molecule, which is only true
    # on a CONSENSUS bam. Here duplicates are merely FLAGGED, not collapsed, so
    # molecule counts are inflated. Use -u for UMI assays, or extend
    # interrogate.py to skip BAM_FDUP reads for non-consensus input.
  fi
  samtools index "$bam"

  echo "[$srr] interrogate"
  python3 "$BINDIR/interrogate.py" run --bam "$bam" --panel "$PANEL" --out "$out"

  [[ $KEEP -eq 1 ]] || rm -rf "$wd"
  echo "[$srr] done -> $out"
}

# simple background-job pool (bash 4.3+: wait -n)
pool() { while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done; }

# --- drive over the manifest (skip header, column 1 = run accession) ---
mapfile -t RUNS < <(tail -n +2 "$MANIFEST" | cut -f1 | awk 'NF')
echo "Manifest: ${#RUNS[@]} runs | jobs=$JOBS threads=$THREADS umi=$UMI -> $OUTDIR"

for srr in "${RUNS[@]}"; do
  if [[ $JOBS -gt 1 ]]; then
    ( process_one "$srr" || echo "$srr" >> "$FAILLOG" ) &
    pool
  else
    process_one "$srr" || echo "$srr" >> "$FAILLOG"
  fi
done
wait

# --- summary ---
done_n=$(find "$OUTDIR" -name '*.tsv' | wc -l | tr -d ' ')
fail_n=$(wc -l < "$FAILLOG" | tr -d ' ')
echo "----"
echo "Outputs: $done_n / ${#RUNS[@]} runs in $OUTDIR"
if [[ $fail_n -gt 0 ]]; then
  echo "FAILED ($fail_n): see $FAILLOG"; exit 1
fi
echo "All runs interrogated. Next: build background.tsv from the Bf blanks, then:"
echo "  python3 $BINDIR/validate.py run-real --manifest $MANIFEST \\"
echo "      --background background.tsv --out results/seqc2_validation.json"
