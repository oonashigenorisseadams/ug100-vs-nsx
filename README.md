# UG 100 vs NovaSeq X Whole-Genome Benchmarking

Reproducibility code for: **"Whole-genome benchmarking reveals context-specific error rates in the Ultima UG100 and Illumina NovaSeqX Platforms."**

## Quick Start

```bash
# 1. Set up environment
conda env create -f environment.yml
conda activate ug100-vs-nsx

# 2. Edit config/paths.yml to point to your local data

# 3. Run a figure (e.g., Figure 1)
sbatch scripts/figure_1/run_happy.sh
```

## Repository Structure

```
├── config/
│   └── paths.yml                  # Central path configuration (edit this first)
├── environment.yml                # Conda environment
├── src/
│   ├── __init__.py
│   └── benchmark_utils.py         # Shared utilities (parsing, config)
├── scripts/
│   ├── extract_duplex_reads.sh    # Extract duplex reads from ppmSeq CRAMs
│   ├── downsample_ppmseq.sh       # Downsample + merge ppmSeq CRAMs
│   ├── figure_1/                  # Whole-genome FP+FN (Fig. 1)
│   │   └── run_happy.sh
│   ├── figure_2/                  # HCR vs LCR stratification (Fig. 2)
│   │   └── run_happy_stratified.sh
│   ├── figure_3/                  # HP-length INDEL accuracy (Fig. 3)
│   │   └── run_happy_hp.sh
│   ├── figure_4/                  # Base-calling error by position (Fig. 4)
│   │   ├── run_error_pileup.sh
│   │   └── error_pileup.py        #   Standalone pileup tool (pysam)
│   ├── figure_5/                  # Read-length distributions (Fig. 5)
│   │   └── run_samtools_stats.sh  #   Also produces GCD for Fig. 6B
│   ├── figure_6/                  # Reproducibility (Fig. 6A)
│   │   ├── run_reproducibility.sh
│   │   └── compute_reproducibility.py
│   └── clinical/                  # ClinVar + tandem repeat overlap
│       └── run_clinical_impact.sh
└── results/                       # Created by scripts (not tracked in git)
```

## Workflow

Each figure has a SLURM data-generation script:

| Figure | Description | Script |
|--------|-------------|--------|
| 1 | Whole-genome FP+FN | `figure_1/run_happy.sh` |
| 2 | HCR vs LCR errors | `figure_2/run_happy_stratified.sh` |
| 3 | INDEL vs HP length | `figure_3/run_happy_hp.sh` |
| 4 | Error by read position | `figure_4/run_error_pileup.sh` |
| 5 | Read-length + GC stats | `figure_5/run_samtools_stats.sh` |
| 6A | Reproducibility | `figure_6/run_reproducibility.sh` |
| Clinical | ClinVar + TR overlap | `clinical/run_clinical_impact.sh` |

Figures 1, 2, 3 can run in parallel (independent hap.py calls).
Figure 5 produces data shared with Figure 6B (samtools stats).

All scripts process four dataset types:
- **NovaSeq X 25B** (3 replicates, DRAGEN 4.4)
- **UG 100 standard WGS** (5 replicates, DeepVariant)
- **UG 100 ppmSeq** (2 runs, DeepVariant)
- **UG 100 ppmSeq duplex** (duplex reads extracted from ppmSeq)

Exceptions:
- **Reproducibility (Fig. 6A)** uses 5 NovaSeq X reps (3 × 25B + 2 × 10B)
  and 5 UG 100 reps (standard WGS only). ppmSeq is excluded per Methods,
  as only 2 sequencing runs were available (insufficient for the analysis).
- **Accuracy benchmarking (Figs 1–5)** uses only NovaSeq X 25B replicates.

### ppmSeq Preprocessing

Before running figure scripts on ppmSeq data:

1. **Extract duplex reads** from raw ppmSeq CRAMs:
   ```bash
   bash scripts/extract_duplex_reads.sh <input.cram> <output.cram> <ref.fasta>
   ```

2. **Downsample and merge** ppmSeq CRAMs to ~36x:
   ```bash
   bash scripts/downsample_ppmseq.sh
   ```

3. **Run variant calling** on merged CRAMs via Ultima AWS Ready2Run DeepVariant v1.0.

## Configuration

All paths are configured in `config/paths.yml`. Edit this once to match your
data layout. Shared utility functions in `src/benchmark_utils.py` read this
config and resolve paths automatically.

## Dependencies

- **hap.py v0.3.15** — via Singularity container (`happy_0.3.15.sif`)
- **DRAGEN v4.4** — for NovaSeq X variant calling
- **Ultima AWS Ready2Run DeepVariant v1.0** — for UG 100 variant calling
- **samtools**, **bcftools**, **bedtools** — see `environment.yml`
- **pysam** — for Figure 4 error pileup

## Data

- Sample: HG002 (NA24385), GRCh38
- Truth: GIAB NIST v4.2.1
- NovaSeq X 25B: 3 replicates, 25B flowcell, DRAGEN 4.4, ~35x
- NovaSeq X 10B: 2 replicates, 10B flowcell, DRAGEN 4.4, ~35x (reproducibility only)
- UG 100: 5 replicates (5 runs), standard WGS chemistry, DeepVariant, ~35x
- UG 100 ppmSeq: 2 runs, ppmSeq chemistry, DeepVariant, ~35x
- UG 100 ppmSeq duplex: duplex reads extracted from ppmSeq (st=et="MIXED"), ~35x

Data access: SRA BioProject PRJNA1427896 [Release 03/13/2026]
