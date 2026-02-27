#!/usr/bin/env bash
# scripts/figure_5/run_samtools_stats.sh
#
# Runs samtools stats on UG CRAMs and NSX BAMs.
# Produces .rl.tsv (read lengths) and .gcd.tsv (GC depth) files
# used by Figure 5 (read length) and Figure 6B (GC bias).
#
#SBATCH --job-name=samtools_stats
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=scripts/logs/samtools_stats_%j.out
#SBATCH --partition=specify

set -euo pipefail

ROOT="EDIT_ME"
OUT_DIR="${ROOT}/results/samtools_stats"
mkdir -p "${OUT_DIR}" scripts/logs

module load samtools 2>/dev/null || true

run_stats() {
    local tag="$1"
    local aln="$2"
    local base
    base=$(basename "${aln}" | sed 's/\.\(cram\|bam\)$//')
    local out_file="${OUT_DIR}/${tag}__${base}.stats"

    if [[ -f "${out_file}" ]]; then
        echo "SKIP: ${base}"
        return
    fi

    echo "Running: ${tag} / ${base}"
    samtools stats \
        -@ 8 \
        --cov-threshold 0 \
        --coverage 0,5000,1 \
        --insert-size 20000 \
        --GC-depth 100 \
        -d -p \
        "${aln}" > "${out_file}"

    grep ^RL  "${out_file}" > "${OUT_DIR}/${tag}__${base}.rl.tsv"  2>/dev/null || true
    grep ^GCD "${out_file}" > "${OUT_DIR}/${tag}__${base}.gcd.tsv" 2>/dev/null || true
    grep ^COV "${out_file}" > "${OUT_DIR}/${tag}__${base}.cov.tsv" 2>/dev/null || true
}

# ── UG merged+downsampled CRAMs (standard WGS) ──────────────
for cram in "${ROOT}"/data/UG_WGS_merged_downsampled/*.cram; do
    [[ -f "${cram}" ]] && run_stats "UG_merged" "${cram}"
done

# ── NovaSeq X BAMs ───────────────────────────────────────────
for bam in "${ROOT}"/data/ILMN_datasets/basespace_novaseqx/bams/NovaSeqX-25B-rep*.bam; do
    [[ -f "${bam}" ]] && run_stats "NSX" "${bam}"
done

# ── UG ppmSeq merged CRAMs ──────────────────────────────────
for cram in "${ROOT}"/data/UG_ppmSeq_merged_downsampled/*.cram; do
    [[ -f "${cram}" ]] && run_stats "UG_ppmSeq" "${cram}"
done

# ── UG ppmSeq duplex merged CRAMs ───────────────────────────
for cram in "${ROOT}"/data/UG_ppmSeq_duplex_merged_downsampled/*.cram; do
    [[ -f "${cram}" ]] && run_stats "UG_ppmSeq_duplex" "${cram}"
done

echo ""
echo "=== Done ==="
