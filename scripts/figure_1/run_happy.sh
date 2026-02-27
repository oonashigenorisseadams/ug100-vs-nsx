#!/usr/bin/env bash
# scripts/figure_1/run_happy.sh
#
# Figure 1: Whole-genome variant calling accuracy (FP+FN).
# Runs hap.py for all callsets against GIAB NIST v4.2.1.
#
# Produces .summary.csv files for downstream analysis.
#
# Usage:
#   sbatch scripts/figure_1/run_happy.sh
#
#SBATCH --job-name=fig1_happy
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=scripts/logs/fig1_happy_%j.out
#SBATCH --partition=specify

set -euo pipefail

# ── Paths (edit config/paths.yml, then source here) ─────────
ROOT="EDIT_ME"

REF_FASTA="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.fna"
REF_SDF="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.sdf"
TRUTH_VCF="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
HAPPY_SIF="${ROOT}/containers/happy_0.3.15.sif"
HAPPY_CMD="singularity exec --bind ${ROOT} ${HAPPY_SIF} /opt/hap.py/bin/hap.py"

THREADS=20

OUT_DIR="${ROOT}/results/fig1/happy"
mkdir -p "${OUT_DIR}" scripts/logs

# ── Callset registry ─────────────────────────────────────────
declare -a CALLSETS=()

# NovaSeq X 25B (DRAGEN 4.4, 35x)
CALLSETS+=(
  "NovaSeqX_rep1|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep1.hard-filtered.vcf.gz"
  "NovaSeqX_rep2|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep2.hard-filtered.vcf.gz"
  "NovaSeqX_rep3|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep3.hard-filtered.vcf.gz"
)

# UG 100 (AWS Ready2Run DeepVariant, ~35x merged)
CALLSETS+=(
  "UG_AWS_DV_411520|${ROOT}/data/callsets/UG_AWS_DeepVariant/411520-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_412157|${ROOT}/data/callsets/UG_AWS_DeepVariant/412157-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_417399|${ROOT}/data/callsets/UG_AWS_DeepVariant/417399-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_418177|${ROOT}/data/callsets/UG_AWS_DeepVariant/418177-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_418399|${ROOT}/data/callsets/UG_AWS_DeepVariant/418399-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
)

# UG 100 ppmSeq (AWS Ready2Run DeepVariant, ~35x merged)
CALLSETS+=(
  "UG_ppmSeq_rep1|${ROOT}/data/callsets/UG_ppmSeq_DeepVariant/ppmSeq_rep1.vcf.gz"
  "UG_ppmSeq_rep2|${ROOT}/data/callsets/UG_ppmSeq_DeepVariant/ppmSeq_rep2.vcf.gz"
)

# UG 100 ppmSeq duplex only (AWS Ready2Run DeepVariant, ~35x merged)
CALLSETS+=(
  "UG_ppmSeq_duplex_rep1|${ROOT}/data/callsets/UG_ppmSeq_duplex_DeepVariant/ppmSeq_duplex_rep1.vcf.gz"
  "UG_ppmSeq_duplex_rep2|${ROOT}/data/callsets/UG_ppmSeq_duplex_DeepVariant/ppmSeq_duplex_rep2.vcf.gz"
)

# ── Sanity checks ────────────────────────────────────────────
for f in "${REF_FASTA}" "${TRUTH_VCF}" "${TRUTH_BED}"; do
    [[ -e "${f}" ]] || { echo "ERROR: Missing ${f}"; exit 1; }
done
[[ -d "${REF_SDF}" ]] || { echo "ERROR: Missing SDF: ${REF_SDF}"; exit 1; }

# ── Index VCFs ───────────────────────────────────────────────
module load bcftools 2>/dev/null || true

for entry in "${CALLSETS[@]}"; do
    label="${entry%%|*}"; vcf="${entry##*|}"
    if [[ ! -f "${vcf}.tbi" ]] && [[ ! -f "${vcf}.csi" ]]; then
        echo "Indexing: ${label}"
        bcftools index -t "${vcf}"
    fi
done

# ── Run hap.py ───────────────────────────────────────────────
for entry in "${CALLSETS[@]}"; do
    label="${entry%%|*}"; vcf="${entry##*|}"
    out_prefix="${OUT_DIR}/${label}_vs_HG002"

    [[ -f "${out_prefix}.summary.csv" ]] && { echo "SKIP: ${label}"; continue; }

    echo "--- ${label} ---"
    ${HAPPY_CMD} \
        "${TRUTH_VCF}" "${vcf}" \
        -r "${REF_FASTA}" -f "${TRUTH_BED}" \
        -o "${out_prefix}" \
        --threads "${THREADS}" --write-counts \
        --engine vcfeval --gender none

    echo "  Done → ${out_prefix}.summary.csv"
done

echo ""
echo "=== Complete ==="
