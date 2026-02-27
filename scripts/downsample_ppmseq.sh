#!/usr/bin/env bash
# scripts/downsample_ppmseq.sh
#
# Downsample and merge ppmSeq CRAMs to ~36x raw coverage.
# Follows the same strategy as standard WGS (Supplementary Note 4):
#   1. Compute raw coverage per library (excl. unmapped/secondary/supplementary)
#   2. Downsample each to 18x
#   3. Merge pairs to reach ~36x
#
# The same approach is used for both all-read ppmSeq and duplex-only ppmSeq.
# For duplex-only: first extract duplex reads with extract_duplex_reads.sh,
# then run this script on the duplex CRAMs.
#
# Usage:
#   bash scripts/downsample_ppmseq.sh
#
# Edit the paths below before running.

set -euo pipefail

# ============================================================
# Configuration — edit these paths
# ============================================================
ROOT="EDIT_ME"
REF_FASTA="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.fna"
AUTOSOMES_BED="${ROOT}/data/ref/GRCh38_autosomes_gaps_removed.bed"

# Input: raw ppmSeq CRAMs (all reads)
PPMSEQ_INPUT_DIR="${ROOT}/data/UG_ppmSeq_raw"
PPMSEQ_OUTPUT_DIR="${ROOT}/data/UG_ppmSeq_merged_downsampled"

# Input: duplex-only ppmSeq CRAMs (from extract_duplex_reads.sh)
DUPLEX_INPUT_DIR="${ROOT}/data/UG_ppmSeq_duplex_raw"
DUPLEX_OUTPUT_DIR="${ROOT}/data/UG_ppmSeq_duplex_merged_downsampled"

TARGET_HALF_COV=18  # Each library downsampled to 18x, pairs merged to ~36x

mkdir -p "${PPMSEQ_OUTPUT_DIR}" "${DUPLEX_OUTPUT_DIR}"

module load samtools 2>/dev/null || true

# ============================================================
# Functions
# ============================================================

compute_mean_coverage() {
    local cram="$1"
    local cov_file="/tmp/$(basename ${cram} .cram)_cov.tsv"

    samtools stats -@ 8 \
        --cov-threshold 0 --coverage 0,5000,1 \
        --filtering-flag 2820 \
        --target-regions "${AUTOSOMES_BED}" \
        "${cram}" | grep ^COV | cut -f 2- > "${cov_file}"

    python3 -c "
import csv
total_bases = total_genome = 0
with open('${cov_file}') as f:
    for row in csv.reader(f, delimiter='\t'):
        cov, nbases = int(row[1]), int(row[2])
        total_bases += cov * nbases
        total_genome += nbases
print(f'{total_bases/total_genome:.2f}')
"
}

downsample_cram() {
    local cram="$1"
    local output="$2"
    local target_cov="$3"

    local mean_cov
    mean_cov=$(compute_mean_coverage "${cram}")

    local frac
    frac=$(python3 -c "print(f'{${target_cov}/${mean_cov}:.6f}')")

    echo "  $(basename ${cram}): mean_cov=${mean_cov}x, frac=${frac}"
    samtools view -s "${frac}" -C \
        --reference "${REF_FASTA}" \
        -o "${output}" "${cram}"
    samtools index "${output}"
}

# ============================================================
# Process ppmSeq all-reads CRAMs
# ============================================================
echo "=== Processing ppmSeq all-reads CRAMs ==="

for cram in "${PPMSEQ_INPUT_DIR}"/*.cram; do
    [[ -f "${cram}" ]] || continue
    base=$(basename "${cram}" .cram)
    ds_cram="${PPMSEQ_OUTPUT_DIR}/${base}_18x.cram"

    if [[ -f "${ds_cram}" ]]; then
        echo "SKIP: ${base} (already downsampled)"
        continue
    fi

    echo "Downsampling: ${base}"
    downsample_cram "${cram}" "${ds_cram}" "${TARGET_HALF_COV}"
done

echo ""
echo "Downsampled ppmSeq CRAMs are in: ${PPMSEQ_OUTPUT_DIR}"
echo "Merge pairs manually with:"
echo "  samtools merge <cram1_18x> <cram2_18x> -o <merged_36x.cram> -O CRAM"
echo "Then run variant calling via Ultima AWS Ready2Run DeepVariant v1.0."

# ============================================================
# Process ppmSeq duplex CRAMs
# ============================================================
echo ""
echo "=== Processing ppmSeq duplex CRAMs ==="

for cram in "${DUPLEX_INPUT_DIR}"/*.cram; do
    [[ -f "${cram}" ]] || continue
    base=$(basename "${cram}" .cram)
    ds_cram="${DUPLEX_OUTPUT_DIR}/${base}_18x.cram"

    if [[ -f "${ds_cram}" ]]; then
        echo "SKIP: ${base} (already downsampled)"
        continue
    fi

    echo "Downsampling: ${base}"
    downsample_cram "${cram}" "${ds_cram}" "${TARGET_HALF_COV}"
done

echo ""
echo "Downsampled duplex CRAMs are in: ${DUPLEX_OUTPUT_DIR}"
echo "Merge pairs and call variants as above."

echo ""
echo "=== Done ==="
