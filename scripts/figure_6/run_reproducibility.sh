#!/usr/bin/env bash
# scripts/figure_6/run_reproducibility.sh
#
# Figure 6A: Variant calling reproducibility across replicates
#
# Per Supplementary Note 6:
#   - Extract PASS SNP and INDEL calls separately
#   - Use bedtools multiinter to find overlaps across 2-5 reps
#   - Compute % reproducible for all combinations of k reps
#   - Repeat genome-wide and restricted to NIST v4.2.1
#
# NovaSeq X uses 5 reps: 3 x 25B + 2 x 10B (Supplementary Note 6)
# UG 100 uses 5 reps: standard WGS chemistry only
#
# NOTE: ppmSeq excluded per Methods — only 2 sequencing runs were
# available, insufficient for meaningful reproducibility analysis
# (the analysis requires ≥3 replicates from independent runs).
#
# Requires: bcftools, bedtools
#
# Usage:
#   sbatch scripts/figure_6/run_reproducibility.sh
#
#SBATCH --job-name=fig6a_repro
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=scripts/logs/fig6a_repro_%j.out
#SBATCH --partition=specify

set -euo pipefail

ROOT="EDIT_ME"
OUT_DIR="${ROOT}/results/fig6a"
mkdir -p "${OUT_DIR}/vcfs" "${OUT_DIR}/multiinter" scripts/logs

TRUTH_BED="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

module load bcftools 2>/dev/null || true
module load bedtools 2>/dev/null || true

# ============================================================
# VCF registry
# ============================================================

# UG 100 — 5 reps (standard WGS chemistry only; ppmSeq excluded per Methods)
declare -a UG_VCFS=(
  "${ROOT}/data/callsets/UG_AWS_DeepVariant/411520-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "${ROOT}/data/callsets/UG_AWS_DeepVariant/412157-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "${ROOT}/data/callsets/UG_AWS_DeepVariant/417399-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "${ROOT}/data/callsets/UG_AWS_DeepVariant/418177-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "${ROOT}/data/callsets/UG_AWS_DeepVariant/418399-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
)

# NovaSeq X — 3 x 25B + 2 x 10B = 5 reps total (Supplementary Note 6)
declare -a NSX_VCFS=(
  "${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep1.hard-filtered.vcf.gz"
  "${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep2.hard-filtered.vcf.gz"
  "${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep3.hard-filtered.vcf.gz"
  "${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-10B-rep1.hard-filtered.vcf.gz"
  "${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-10B-rep2.hard-filtered.vcf.gz"
)

# ============================================================
# STEP 1: Extract PASS SNPs and INDELs from each VCF
# ============================================================
echo "=== Extracting PASS SNPs and INDELs ==="

extract_pass() {
    local vcf="$1"
    local label="$2"
    local vtype="$3"  # snps or indels

    local out="${OUT_DIR}/vcfs/${label}_${vtype}.vcf.gz"
    if [[ -f "${out}" ]]; then
        echo "  OK: ${label} ${vtype}"
        return
    fi

    echo "  Extracting: ${label} ${vtype}"
    bcftools view -f PASS -v "${vtype}" "${vcf}" -Oz -o "${out}"
    bcftools index -t "${out}"
}

# UG
for i in "${!UG_VCFS[@]}"; do
    vcf="${UG_VCFS[$i]}"
    run_id=$(basename "${vcf}" | cut -d'-' -f1)
    extract_pass "${vcf}" "UG_${run_id}" "snps"
    extract_pass "${vcf}" "UG_${run_id}" "indels"
done

# NSX
for i in "${!NSX_VCFS[@]}"; do
    vcf="${NSX_VCFS[$i]}"
    rep=$((i + 1))
    extract_pass "${vcf}" "NSX_rep${rep}" "snps"
    extract_pass "${vcf}" "NSX_rep${rep}" "indels"
done

# ============================================================
# STEP 2: Run bedtools multiinter for each platform x variant type
# ============================================================
echo ""
echo "=== Running bedtools multiinter ==="

run_multiinter() {
    local platform="$1"
    local vtype="$2"
    shift 2
    local vcfs=("$@")

    local out="${OUT_DIR}/multiinter/${platform}_${vtype}_multiinter.bed"
    local out_nist="${OUT_DIR}/multiinter/${platform}_${vtype}_multiinter_NIST.bed"

    if [[ -f "${out}" ]] && [[ -f "${out_nist}" ]]; then
        echo "  OK: ${platform} ${vtype}"
        return
    fi

    echo "  ${platform} ${vtype}: ${#vcfs[@]} VCFs"
    bedtools multiinter -header -i "${vcfs[@]}" > "${out}"

    # Restrict to NIST regions
    bedtools intersect -a "${out}" -b "${TRUTH_BED}" -header > "${out_nist}"
}

# Collect extracted VCFs for each platform x type
UG_SNP_VCFS=()
UG_INDEL_VCFS=()
for i in "${!UG_VCFS[@]}"; do
    run_id=$(basename "${UG_VCFS[$i]}" | cut -d'-' -f1)
    UG_SNP_VCFS+=("${OUT_DIR}/vcfs/UG_${run_id}_snps.vcf.gz")
    UG_INDEL_VCFS+=("${OUT_DIR}/vcfs/UG_${run_id}_indels.vcf.gz")
done

NSX_SNP_VCFS=()
NSX_INDEL_VCFS=()
for i in "${!NSX_VCFS[@]}"; do
    rep=$((i + 1))
    NSX_SNP_VCFS+=("${OUT_DIR}/vcfs/NSX_rep${rep}_snps.vcf.gz")
    NSX_INDEL_VCFS+=("${OUT_DIR}/vcfs/NSX_rep${rep}_indels.vcf.gz")
done

run_multiinter "UG" "snps" "${UG_SNP_VCFS[@]}"
run_multiinter "UG" "indels" "${UG_INDEL_VCFS[@]}"
run_multiinter "NSX" "snps" "${NSX_SNP_VCFS[@]}"
run_multiinter "NSX" "indels" "${NSX_INDEL_VCFS[@]}"

# ============================================================
# STEP 3: Compute reproducibility stats
# ============================================================
echo ""
echo "=== Computing reproducibility ==="

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
python3 "${SCRIPT_DIR}/compute_reproducibility.py" \
    --multiinter-dir "${OUT_DIR}/multiinter" \
    --output "${OUT_DIR}/fig6a_reproducibility.csv"

echo ""
echo "=== Done ==="
echo "Results in: ${OUT_DIR}/"
