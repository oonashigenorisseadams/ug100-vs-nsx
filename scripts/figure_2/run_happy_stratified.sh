#!/usr/bin/env bash
# scripts/figure_2/run_happy_stratified.sh
#
# Figure 2: FP+FN stratified by UG HCR vs UG LCR.
# Runs hap.py with --stratification for UG HCR/LCR regions.
#
#SBATCH --job-name=fig2_happy
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=scripts/logs/fig2_happy_%j.out
#SBATCH --partition=specify

set -euo pipefail

ROOT="EDIT_ME"

REF_FASTA="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.fna"
TRUTH_VCF="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
UG_HCR_BED="${ROOT}/data/ref/UG_HCR_v3.1/ug_hcr.bed"

HAPPY_CMD="singularity exec --bind ${ROOT} ${ROOT}/containers/happy_0.3.15.sif /opt/hap.py/bin/hap.py"
THREADS=20

STRAT_DIR="${ROOT}/results/fig2/stratification"
OUT_DIR="${ROOT}/results/fig2/happy"
mkdir -p "${STRAT_DIR}" "${OUT_DIR}" scripts/logs

# ── Prepare HCR/LCR stratification beds ──────────────────────
module load bedtools 2>/dev/null || true
module load bcftools 2>/dev/null || true

UG_HCR_SORTED="${STRAT_DIR}/ug_hcr_merged_sorted.bed"
UG_LCR_BED="${STRAT_DIR}/ug_lcr.bed"
STRAT_TSV="${STRAT_DIR}/hcr_stratification.tsv"

if [[ ! -s "${UG_HCR_SORTED}" ]]; then
    bedtools sort -i "${UG_HCR_BED}" -g "${REF_FASTA}.fai" | bedtools merge -i stdin > "${UG_HCR_SORTED}"
fi

if [[ ! -f "${UG_LCR_BED}" ]]; then
    bedtools complement -i "${UG_HCR_SORTED}" -g "${REF_FASTA}.fai" > "${UG_LCR_BED}"
fi

echo -e "UG_HCR\t${UG_HCR_SORTED}" > "${STRAT_TSV}"
echo -e "UG_LCR\t${UG_LCR_BED}" >> "${STRAT_TSV}"

# ── Callset registry ─────────────────────────────────────────
declare -a CALLSETS=()

# NovaSeq X 25B
CALLSETS+=(
  "NovaSeqX_rep1|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep1.hard-filtered.vcf.gz"
  "NovaSeqX_rep2|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep2.hard-filtered.vcf.gz"
  "NovaSeqX_rep3|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep3.hard-filtered.vcf.gz"
)

# UG 100 standard WGS
CALLSETS+=(
  "UG_AWS_DV_411520|${ROOT}/data/callsets/UG_AWS_DeepVariant/411520-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_412157|${ROOT}/data/callsets/UG_AWS_DeepVariant/412157-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_417399|${ROOT}/data/callsets/UG_AWS_DeepVariant/417399-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_418177|${ROOT}/data/callsets/UG_AWS_DeepVariant/418177-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
  "UG_AWS_DV_418399|${ROOT}/data/callsets/UG_AWS_DeepVariant/418399-RM8391-HG002-1-Z0105-B-CATGGCGCAGTGCTGAT_Z0106-B-CTGAGAATGTGTGAT.cram.vcf.gz"
)

# UG 100 ppmSeq
CALLSETS+=(
  "UG_ppmSeq_rep1|${ROOT}/data/callsets/UG_ppmSeq_DeepVariant/ppmSeq_rep1.vcf.gz"
  "UG_ppmSeq_rep2|${ROOT}/data/callsets/UG_ppmSeq_DeepVariant/ppmSeq_rep2.vcf.gz"
)

# UG 100 ppmSeq duplex only
CALLSETS+=(
  "UG_ppmSeq_duplex_rep1|${ROOT}/data/callsets/UG_ppmSeq_duplex_DeepVariant/ppmSeq_duplex_rep1.vcf.gz"
  "UG_ppmSeq_duplex_rep2|${ROOT}/data/callsets/UG_ppmSeq_duplex_DeepVariant/ppmSeq_duplex_rep2.vcf.gz"
)

# ── Index + Run ──────────────────────────────────────────────
for entry in "${CALLSETS[@]}"; do
    label="${entry%%|*}"; vcf="${entry##*|}"
    [[ -f "${vcf}.tbi" ]] || [[ -f "${vcf}.csi" ]] || bcftools index -t "${vcf}"
done

for entry in "${CALLSETS[@]}"; do
    label="${entry%%|*}"; vcf="${entry##*|}"
    out_prefix="${OUT_DIR}/${label}_vs_HG002_hcr"

    [[ -f "${out_prefix}.extended.csv" ]] && { echo "SKIP: ${label}"; continue; }

    echo "--- ${label} ---"
    ${HAPPY_CMD} \
        "${TRUTH_VCF}" "${vcf}" \
        -r "${REF_FASTA}" -f "${TRUTH_BED}" \
        -o "${out_prefix}" \
        --threads "${THREADS}" --write-counts \
        --engine vcfeval --gender none \
        --stratification "${STRAT_TSV}"
done

echo "=== Complete ==="
