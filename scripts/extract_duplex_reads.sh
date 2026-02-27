#!/usr/bin/env bash
# scripts/extract_duplex_reads.sh
#
# Extract duplex reads from ppmSeq CRAMs.
# Duplex reads are identified by st and et tags both set to "MIXED"
# (per manuscript Methods).
#
# Usage:
#   bash scripts/extract_duplex_reads.sh <input.cram> <output.cram> <ref.fasta>
#
# Example:
#   bash scripts/extract_duplex_reads.sh \
#       data/UG_ppmSeq_raw/run1_rep1.cram \
#       data/UG_ppmSeq_duplex_raw/run1_rep1_duplex.cram \
#       data/ref/GRCh38_no_alt_analysis_set.fna

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input.cram> <output.cram> <ref.fasta>"
    exit 1
fi

INPUT_CRAM="$1"
OUTPUT_CRAM="$2"
REF_FASTA="$3"

[[ -f "${INPUT_CRAM}" ]] || { echo "ERROR: Input not found: ${INPUT_CRAM}"; exit 1; }
[[ -f "${REF_FASTA}" ]]  || { echo "ERROR: Reference not found: ${REF_FASTA}"; exit 1; }

mkdir -p "$(dirname "${OUTPUT_CRAM}")"

echo "Extracting duplex reads (st=MIXED, et=MIXED) from $(basename ${INPUT_CRAM})..."
samtools view \
    -@ 8 \
    -C \
    --reference "${REF_FASTA}" \
    -e '[st] == "MIXED" && [et] == "MIXED"' \
    -o "${OUTPUT_CRAM}" \
    "${INPUT_CRAM}"

echo "Indexing output..."
samtools index "${OUTPUT_CRAM}"

echo "Done: ${OUTPUT_CRAM}"
echo "  Input reads:  $(samtools view -c "${INPUT_CRAM}")"
echo "  Duplex reads: $(samtools view -c "${OUTPUT_CRAM}")"
