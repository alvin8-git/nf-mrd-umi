#!/usr/bin/env bash
# build_panel.sh
#
# Build the tumor-informed panel for the SEQC2 cell-line validation cohort:
#
#   SEQC2 high-confidence HCC1395 somatic SNVs
#     INTERSECT  high-confidence regions
#     INTERSECT  the assay's target BED
#   -> panel.vcf.gz (+ panel.bed)
#
# In a real patient, Pipeline A (panel_select.py) produces this from the
# patient's tumor WES. For the SEQC2 HCC1395/HCC1395BL cohort the "tumor"
# mutations are already known (the SEQC2 truth set), so the panel is just the
# truth SNVs the assay can actually see. Output feeds interrogate.py --panel.
#
# Truth set: SEQC2 Somatic Mutation WG (Fang et al. 2021), BioProject
# PRJNA489865, files verified at the NCBI reference-samples FTP below.

set -euo pipefail

TRUTH_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest"
VCF_URL="$TRUTH_BASE/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
HCBED_URL="$TRUTH_BASE/High-Confidence_Regions_v1.2.bed"

usage() {
  cat <<EOF
Usage: $(basename "$0") -b ASSAY_TARGET.bed [options]

Required:
  -b FILE   assay target BED (the regions the ctDNA assay covers, e.g. TFS2).
            Source: SEQC2 Liquid Biopsy figshare (panel region BEDs per vendor).

Options:
  -o PREFIX output prefix              (default: panel)  -> PREFIX.vcf.gz, PREFIX.bed
  -v FILE   SEQC2 high-conf SNV VCF    (default: download verified NCBI file)
  -c FILE   high-confidence regions BED(default: download verified NCBI file)
  -r FILE   reference FASTA            (optional: left-normalize + chrom check)
  -d DIR    download cache dir         (default: ref/seqc2)
  -C STYLE  rename chroms to match ref: 'ucsc' (chr1) or 'ensembl' (1) (default: none)
  -h        this help

Truth set defaults (verified to exist):
  VCF : $VCF_URL
  BED : $HCBED_URL
EOF
}

PREFIX="panel"; ASSAY_BED=""; SNV_VCF=""; HC_BED=""; REF=""; DLDIR="ref/seqc2"; CHRSTYLE=""
while getopts "b:o:v:c:r:d:C:h" opt; do
  case "$opt" in
    b) ASSAY_BED=$OPTARG ;; o) PREFIX=$OPTARG ;; v) SNV_VCF=$OPTARG ;;
    c) HC_BED=$OPTARG ;; r) REF=$OPTARG ;; d) DLDIR=$OPTARG ;;
    C) CHRSTYLE=$OPTARG ;; h) usage; exit 0 ;; *) usage; exit 2 ;;
  esac
done

[[ -n $ASSAY_BED && -f $ASSAY_BED ]] || { echo "ERROR: -b assay target BED missing" >&2; usage; exit 2; }
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing tool: $1" >&2; exit 2; }; }
need bcftools; need bedtools; need tabix; need bgzip

fetch() { # url dest
  [[ -s $2 ]] && { echo "  cached: $2"; return 0; }
  if command -v curl >/dev/null 2>&1; then curl -fSL "$1" -o "$2"
  elif command -v wget >/dev/null 2>&1; then wget -O "$2" "$1"
  else echo "ERROR: need curl or wget to download truth set" >&2; exit 2; fi
}

mkdir -p "$DLDIR"
if [[ -z $SNV_VCF ]]; then SNV_VCF="$DLDIR/$(basename "$VCF_URL")"; echo "Fetch truth SNVs:"; fetch "$VCF_URL" "$SNV_VCF"; fi
if [[ -z $HC_BED ]];  then HC_BED="$DLDIR/$(basename "$HCBED_URL")"; echo "Fetch HC regions:"; fetch "$HCBED_URL" "$HC_BED"; fi
[[ -f ${SNV_VCF}.tbi ]] || tabix -p vcf "$SNV_VCF" 2>/dev/null || true

# Optional chrom-name harmonization so VCF/BED match the alignment reference.
# ponytail: simple sed on the BEDs + bcftools annotate on the VCF; covers the
# common chr<->no-chr mismatch, not exotic contig naming.
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
A_BED="$ASSAY_BED"; C_BED="$HC_BED"; V_VCF="$SNV_VCF"
if [[ $CHRSTYLE == ucsc || $CHRSTYLE == ensembl ]]; then
  echo "Harmonizing chrom names to: $CHRSTYLE"
  if [[ $CHRSTYLE == ucsc ]]; then
    addchr() { sed -E 's/^([0-9XYM])/chr\1/' "$1"; }
  else
    addchr() { sed -E 's/^chr//' "$1"; }
  fi
  A_BED="$work/assay.bed"; C_BED="$work/hc.bed"
  addchr "$ASSAY_BED" > "$A_BED"; addchr "$HC_BED" > "$C_BED"
  # rename VCF contigs via a map file
  bcftools view -h "$SNV_VCF" | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*/\1/' \
    | while read -r c; do
        if [[ $CHRSTYLE == ucsc ]]; then n=${c#chr}; echo "$c chr${n}"; else echo "$c ${c#chr}"; fi
      done > "$work/chrmap.txt"
  V_VCF="$work/snv.renamed.vcf.gz"
  bcftools annotate --rename-chrs "$work/chrmap.txt" -Oz -o "$V_VCF" "$SNV_VCF"
  tabix -p vcf "$V_VCF"
fi

# 1) panel regions = assay BED INTERSECT high-confidence regions
echo "Intersecting assay BED with high-confidence regions..."
bedtools intersect -a "$A_BED" -b "$C_BED" \
  | sort -k1,1 -k2,2n | bedtools merge > "$work/panel_regions.bed"
nreg=$(wc -l < "$work/panel_regions.bed" | tr -d ' ')
echo "  panel regions: $nreg intervals"

# 2) select PASS, biallelic SNVs within those regions; optionally left-normalize
echo "Selecting biallelic PASS SNVs in panel regions..."
norm=(cat)
[[ -n $REF ]] && norm=(bcftools norm -f "$REF" -Ou -)
bcftools view -R "$work/panel_regions.bed" -v snps -m2 -M2 -f PASS,. -Ou "$V_VCF" \
  | "${norm[@]}" \
  | bcftools sort -Oz -o "${PREFIX}.vcf.gz"
tabix -p vcf "${PREFIX}.vcf.gz"

# 3) companion BED (0-based) for interrogate / QC
bcftools query -f '%CHROM\t%POS0\t%END\t%REF>%ALT\n' "${PREFIX}.vcf.gz" > "${PREFIX}.bed"

n=$(bcftools view -H "${PREFIX}.vcf.gz" | wc -l | tr -d ' ')
echo "----"
echo "Panel: $n SNV sites -> ${PREFIX}.vcf.gz (+ ${PREFIX}.bed)"
if [[ $n -lt 8 ]]; then
  echo "WARN: only $n trackable sites. Panel-integrated MRD wants more (see TODO.md"
  echo "      'minimum-trackable-N gate'). A small hotspot assay may simply not"
  echo "      cover enough HCC1395 mutations; consider a broader assay's BED."
fi
echo "Next: bin/fetch_and_interrogate.sh -p ${PREFIX}.vcf.gz ..."
