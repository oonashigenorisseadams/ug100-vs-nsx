#!/usr/bin/env bash
# scripts/figure_4/run_error_pileup.sh
#
# Figure 4: Base-calling error rates by position in read.
# Runs fig4_error_pileup.py on each BAM/CRAM.
#
# Per Supplementary Note 3:
#   - MAPQ > 30, within NIST v4.2.1 high-confidence regions
#   - Mask truth VCF variants
#   - Count mismatches, insertions, deletions by read position
#
#SBATCH --job-name=fig4_err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=scripts/logs/fig4_err_%j.out
#SBATCH --partition=specify

set -euo pipefail

ROOT="EDIT_ME"
OUT_DIR="${ROOT}/results/fig4"
mkdir -p "${OUT_DIR}" scripts/logs

REF_FASTA="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.fna"
TRUTH_VCF="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

NSX_BAM_DIR="${ROOT}/data/ILMN_datasets/basespace_novaseqx/bams"
UG_CRAM_DIR="${ROOT}/data/UG_WGS_merged_downsampled"
UG_PPMSEQ_CRAM_DIR="${ROOT}/data/UG_ppmSeq_merged_downsampled"
UG_PPMSEQ_DUPLEX_CRAM_DIR="${ROOT}/data/UG_ppmSeq_duplex_merged_downsampled"

# ── Build truth variant mask ─────────────────────────────────
TRUTH_MASK="${OUT_DIR}/truth_variant_positions.bed.gz"
if [[ ! -f "${TRUTH_MASK}" ]]; then
    module load bcftools 2>/dev/null || true
    module load bedtools 2>/dev/null || true
    bcftools query -f '%CHROM\t%POS0\t%END\n' "${TRUTH_VCF}" \
        | bedtools intersect -a stdin -b "${TRUTH_BED}" \
        | bedtools sort -i stdin -g "${REF_FASTA}.fai" \
        | bedtools merge -i stdin \
        | gzip > "${TRUTH_MASK}"
fi

# ── Run pileup on each alignment ────────────────────────────
PILEUP_SCRIPT="$(dirname "$0")/error_pileup.py"

# NovaSeq X BAMs
for bam in "${NSX_BAM_DIR}"/NovaSeqX-25B-rep*.bam; do
    [[ -f "${bam}" ]] || continue
    label=$(basename "${bam}" .bam)
    out_tsv="${OUT_DIR}/${label}_error_by_pos.tsv.gz"
    [[ -f "${out_tsv}" ]] && { echo "SKIP: ${label}"; continue; }

    echo "--- ${label} ---"
    python3 "${PILEUP_SCRIPT}" \
        --bam "${bam}" --ref "${REF_FASTA}" \
        --regions "${TRUTH_BED}" --mask "${TRUTH_MASK}" \
        --output "${out_tsv}" \
        --platform "NovaSeqX_25B" --label "${label}" \
        --max-position 350
done

# UG 100 standard WGS CRAMs
for cram in "${UG_CRAM_DIR}"/*.cram; do
    [[ -f "${cram}" ]] || continue
    run_id=$(basename "${cram}" | cut -d'-' -f1)
    label="UG_100_${run_id}"
    out_tsv="${OUT_DIR}/${label}_error_by_pos.tsv.gz"
    [[ -f "${out_tsv}" ]] && { echo "SKIP: ${label}"; continue; }

    echo "--- ${label} ---"
    python3 "${PILEUP_SCRIPT}" \
        --bam "${cram}" --ref "${REF_FASTA}" \
        --regions "${TRUTH_BED}" --mask "${TRUTH_MASK}" \
        --output "${out_tsv}" \
        --platform "UG_100" --label "${label}" \
        --max-position 350
done

# UG 100 ppmSeq CRAMs
for cram in "${UG_PPMSEQ_CRAM_DIR}"/*.cram; do
    [[ -f "${cram}" ]] || continue
    run_id=$(basename "${cram}" | cut -d'-' -f1)
    label="UG_100_ppmSeq_${run_id}"
    out_tsv="${OUT_DIR}/${label}_error_by_pos.tsv.gz"
    [[ -f "${out_tsv}" ]] && { echo "SKIP: ${label}"; continue; }

    echo "--- ${label} ---"
    python3 "${PILEUP_SCRIPT}" \
        --bam "${cram}" --ref "${REF_FASTA}" \
        --regions "${TRUTH_BED}" --mask "${TRUTH_MASK}" \
        --output "${out_tsv}" \
        --platform "UG_100_ppmSeq" --label "${label}" \
        --max-position 350
done

# UG 100 ppmSeq duplex CRAMs
for cram in "${UG_PPMSEQ_DUPLEX_CRAM_DIR}"/*.cram; do
    [[ -f "${cram}" ]] || continue
    run_id=$(basename "${cram}" | cut -d'-' -f1)
    label="UG_100_ppmSeq_duplex_${run_id}"
    out_tsv="${OUT_DIR}/${label}_error_by_pos.tsv.gz"
    [[ -f "${out_tsv}" ]] && { echo "SKIP: ${label}"; continue; }

    echo "--- ${label} ---"
    python3 "${PILEUP_SCRIPT}" \
        --bam "${cram}" --ref "${REF_FASTA}" \
        --regions "${TRUTH_BED}" --mask "${TRUTH_MASK}" \
        --output "${out_tsv}" \
        --platform "UG_100_ppmSeq_duplex" --label "${label}" \
        --max-position 350
done

echo "=== Complete ==="
