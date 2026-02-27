#!/usr/bin/env bash
# scripts/clinical/run_clinical_impact.sh
#
# Clinical impact assessment (Supplementary Note 8):
#   1. % ClinVar P+LP (≥2 star) variant sites outside UG HCR v3.1
#   2. % Tandem repeat catalog bases outside UG HCR v3.1
#
# Fast — just bedtools intersections, no heavy compute.
# Can run immediately, no dependencies on other figures.
#
# Usage:
#   sbatch scripts/clinical/run_clinical_impact.sh
#
#SBATCH --job-name=clinical
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=scripts/logs/clinical_%j.out
#SBATCH --partition=specify

set -euo pipefail

# ============================================================
# Paths — edit config/paths.yml or override here
# ============================================================
ROOT="EDIT_ME"

REF_FAI="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.fna.fai"
UG_HCR_BED="${ROOT}/data/ref/UG_HCR_v3.1/ug_hcr.bed"
CLINVAR_VCF="${ROOT}/data/clinvar/clinvar_20240805.vcf.gz"
TR_BED="${ROOT}/data/tandem_repeat_catalog/repeat_catalog_v1.hg38.1_to_1000bp_motifs.TRGT.bed"

OUT_DIR="${ROOT}/results/clinical"
mkdir -p "${OUT_DIR}" scripts/logs

module load bcftools 2>/dev/null || true
module load bedtools 2>/dev/null || true

# ============================================================
# STEP 0: Prepare sorted/merged UG HCR
# ============================================================
echo "=== Preparing UG HCR ==="
UG_HCR_SORTED="${OUT_DIR}/ug_hcr_merged_sorted.bed"

if [[ ! -f "${UG_HCR_SORTED}" ]]; then
    cat "${UG_HCR_BED}" \
        | bedtools sort -i stdin \
        | bedtools merge -i stdin \
        > "${UG_HCR_SORTED}"
    echo "  Wrote ${UG_HCR_SORTED}"
fi

HCR_BP=$(awk '{SUM+=$3-$2} END {print SUM}' "${UG_HCR_SORTED}")
GENOME_BP=$(awk '{SUM+=$2} END {print SUM}' "${REF_FAI}")
HCR_PCT=$(echo "scale=2; ${HCR_BP} * 100 / ${GENOME_BP}" | bc)
echo "  UG HCR: ${HCR_BP} bp = ${HCR_PCT}% of genome (expect ~90.3%)"

# ============================================================
# STEP 1: ClinVar P+LP ≥2 stars
# ============================================================
echo ""
echo "=== ClinVar P+LP analysis ==="

CLINVAR_FILTERED="${OUT_DIR}/clinvar_P_LP_2star.vcf"
CLINVAR_BED="${OUT_DIR}/clinvar_P_LP_2star_merged_sorted.bed"
CLINVAR_IN_HCR="${OUT_DIR}/clinvar_P_LP_2star_in_ug_hcr.bed"

# Extract P+LP ≥2 star variants
if [[ ! -f "${CLINVAR_FILTERED}" ]]; then
    echo "  Extracting P+LP ≥2 star variants..."
    bcftools view "${CLINVAR_VCF}" \
        -i '((CLNSIG~"Likely_pathogenic") | (CLNSIG~"Pathogenic")) & ((CLNREVSTAT="criteria_provided" & CLNREVSTAT="_multiple_submitters" & CLNREVSTAT="_no_conflicts") | CLNREVSTAT="reviewed_by_expert_panel" | CLNREVSTAT="practice_guideline")' \
        > "${CLINVAR_FILTERED}"
    N_VARS=$(grep -c -v "^#" "${CLINVAR_FILTERED}")
    echo "  Extracted ${N_VARS} P+LP variants"
fi

# Convert to BED, sort, merge
if [[ ! -f "${CLINVAR_BED}" ]]; then
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "${CLINVAR_FILTERED}" \
        | bedtools sort -i stdin \
        | bedtools merge -i stdin \
        > "${CLINVAR_BED}"
    echo "  ClinVar bed: $(wc -l < "${CLINVAR_BED}") intervals"
fi

# Intersect with UG HCR
if [[ ! -f "${CLINVAR_IN_HCR}" ]]; then
    bedtools intersect -a "${UG_HCR_SORTED}" -b "${CLINVAR_BED}" > "${CLINVAR_IN_HCR}"
fi

# Calculate % outside HCR
OVERLAP_BP=$(awk '{SUM+=$3-$2} END {print SUM}' "${CLINVAR_IN_HCR}")
CLINVAR_BP=$(awk '{SUM+=$3-$2} END {print SUM}' "${CLINVAR_BED}")
CLINVAR_PCT=$(echo "scale=2; 100 - ${OVERLAP_BP} * 100 / ${CLINVAR_BP}" | bc)

echo ""
echo "  ClinVar P+LP total bp:          ${CLINVAR_BP}"
echo "  ClinVar P+LP in UG HCR:         ${OVERLAP_BP} bp"
echo "  % ClinVar P+LP NOT in UG HCR:   ${CLINVAR_PCT}%  (expect 2.24%)"

# ============================================================
# STEP 2: Tandem Repeat Catalog
# ============================================================
echo ""
echo "=== Tandem Repeat Catalog analysis ==="

TR_SORTED="${OUT_DIR}/repeat_catalog_merged_sorted.bed"
TR_IN_HCR="${OUT_DIR}/repeat_catalog_in_ug_hcr.bed"

# Sort and merge TR catalog
if [[ ! -f "${TR_SORTED}" ]]; then
    echo "  Sorting/merging tandem repeat catalog..."
    cat "${TR_BED}" \
        | bedtools sort -i stdin \
        | bedtools merge -i stdin \
        > "${TR_SORTED}"
    echo "  TR bed: $(wc -l < "${TR_SORTED}") intervals"
fi

# Intersect with UG HCR
if [[ ! -f "${TR_IN_HCR}" ]]; then
    bedtools intersect -a "${UG_HCR_SORTED}" -b "${TR_SORTED}" > "${TR_IN_HCR}"
fi

# Calculate % outside HCR
TR_OVERLAP_BP=$(awk '{SUM+=$3-$2} END {print SUM}' "${TR_IN_HCR}")
TR_TOTAL_BP=$(awk '{SUM+=$3-$2} END {print SUM}' "${TR_SORTED}")
TR_OUTSIDE_BP=$(echo "${TR_TOTAL_BP} - ${TR_OVERLAP_BP}" | bc)
TR_PCT=$(echo "scale=2; 100 - ${TR_OVERLAP_BP} * 100 / ${TR_TOTAL_BP}" | bc)

echo ""
echo "  TR catalog total bp:           ${TR_TOTAL_BP}  (expect 65,678,112)"
echo "  TR catalog in UG HCR:          ${TR_OVERLAP_BP} bp"
echo "  TR catalog NOT in HCR:         ${TR_OUTSIDE_BP} bp  (expect 14,832,550)"
echo "  % TR catalog NOT in UG HCR:    ${TR_PCT}%  (expect 22.58%)"

# ============================================================
# Summary
# ============================================================
SUMMARY="${OUT_DIR}/clinical_impact_summary.txt"
cat > "${SUMMARY}" <<EOF
Clinical Impact Assessment
==========================
UG HCR v3.1 covers ${HCR_PCT}% of GRCh38

ClinVar P+LP (>=2 star, release 2024-08-05):
  Total sites:      ${CLINVAR_BP} bp
  Outside UG HCR:   ${CLINVAR_PCT}%

Tandem Repeat Catalog:
  Total bases:      ${TR_TOTAL_BP} bp
  Outside UG HCR:   ${TR_OUTSIDE_BP} bp (${TR_PCT}%)
EOF

echo ""
echo "=== Done ==="
echo "Summary: ${SUMMARY}"
