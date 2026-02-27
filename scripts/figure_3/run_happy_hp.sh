#!/usr/bin/env bash
# scripts/figure_3/run_happy_hp.sh
#
# Figure 3: INDEL precision & recall vs homopolymer length.
# Runs hap.py with per-HP-length stratification beds (2-20bp, >20bp).
#
#SBATCH --job-name=fig3_hp
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/logs/fig3_hp_%j.out
#SBATCH --partition=specify

set -euo pipefail

ROOT="EDIT_ME"

REF_FASTA="${ROOT}/data/ref/GRCh38_no_alt_analysis_set.fna"
TRUTH_VCF="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${ROOT}/data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
HP_DIR="${ROOT}/data/HP_beds_ge2/by_len"

HAPPY_CMD="singularity exec --bind ${ROOT} ${ROOT}/containers/happy_0.3.15.sif /opt/hap.py/bin/hap.py"
THREADS=10

OUT_DIR="${ROOT}/results/fig3/happy"
STRAT_DIR="${ROOT}/results/fig3/stratification"
mkdir -p "${OUT_DIR}" "${STRAT_DIR}" scripts/logs

module load bedtools 2>/dev/null || true
module load bcftools 2>/dev/null || true

# ── Merge HP 21-30bp into ">20bp" bed ────────────────────────
GT20_BED="${STRAT_DIR}/HP_gt20bp.bed.gz"
if [[ ! -f "${GT20_BED}" ]]; then
    echo "Merging HP 21-30bp into >20bp bed..."
    GT20_BEDS=""
    for len in $(seq 21 30); do
        padded=$(printf "%04d" ${len})
        bed="${HP_DIR}/HP_len_${padded}.bed.gz"
        [[ -f "${bed}" ]] && GT20_BEDS="${GT20_BEDS} ${bed}"
    done
    zcat ${GT20_BEDS} \
        | bedtools sort -i stdin -g "${REF_FASTA}.fai" \
        | bedtools merge -i stdin \
        | gzip > "${GT20_BED}"
fi

# ── Build stratification TSV (20 strata) ─────────────────────
STRAT_TSV="${STRAT_DIR}/hp_stratification.tsv"
echo -n "" > "${STRAT_TSV}"
for len in $(seq 2 20); do
    padded=$(printf "%04d" ${len})
    bed="${HP_DIR}/HP_len_${padded}.bed.gz"
    [[ -f "${bed}" ]] && echo -e "HP_${len}bp\t${bed}" >> "${STRAT_TSV}"
done
echo -e "HP_gt20bp\t${GT20_BED}" >> "${STRAT_TSV}"
echo "Stratification: $(wc -l < "${STRAT_TSV}") strata"

# ── Callset registry ─────────────────────────────────────────
declare -a CALLSETS=()
CALLSETS+=(
  "NovaSeqX_rep1|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep1.hard-filtered.vcf.gz"
  "NovaSeqX_rep2|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep2.hard-filtered.vcf.gz"
  "NovaSeqX_rep3|${ROOT}/data/ILMN_datasets/basespace_novaseqx/NovaSeqX-25B-rep3.hard-filtered.vcf.gz"
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

# ── Run hap.py ───────────────────────────────────────────────
for entry in "${CALLSETS[@]}"; do
    label="${entry%%|*}"; vcf="${entry##*|}"
    [[ -f "${vcf}.tbi" ]] || [[ -f "${vcf}.csi" ]] || bcftools index -t "${vcf}"
    out_prefix="${OUT_DIR}/${label}_vs_HG002_hp"

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
