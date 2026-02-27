"""
src/benchmark_utils.py — Shared utilities for UG100 vs NovaSeq X benchmarking.

Provides:
  - Config loading from paths.yml
  - hap.py summary.csv and extended.csv parsing
  - Homopolymer stratum parsing
  - samtools stats RL / GCD parsing
  - bedtools multiinter parsing
  - Median computation across replicates
  - Duplex read extraction from ppmSeq CRAMs
  - Consistent color/marker definitions
"""

import csv
import gzip
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Optional

import yaml

# ============================================================
# Config
# ============================================================

_CONFIG_CACHE = None


def load_config(config_path: Optional[Path] = None) -> dict:
    """
    Load paths.yml and resolve all relative paths against project_root.

    Parameters
    ----------
    config_path : Path, optional
        Explicit path to paths.yml.  If None, searches upward from this
        file for config/paths.yml.

    Returns
    -------
    dict with 'root' (Path) plus all other keys from the YAML.
    """
    global _CONFIG_CACHE
    if _CONFIG_CACHE is not None and config_path is None:
        return _CONFIG_CACHE

    if config_path is None:
        # Walk upward from this file to find config/paths.yml
        search = Path(__file__).resolve().parent
        for _ in range(5):
            candidate = search / "config" / "paths.yml"
            if candidate.exists():
                config_path = candidate
                break
            search = search.parent
        if config_path is None:
            raise FileNotFoundError(
                "Could not find config/paths.yml. "
                "Pass config_path explicitly or run from the project root."
            )

    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    root = Path(cfg.get("project_root", ".")).resolve()
    cfg["root"] = root
    _CONFIG_CACHE = cfg
    return cfg


def resolve(cfg: dict, *keys: str) -> Path:
    """
    Resolve a dotted config key to an absolute Path.

    Examples
    --------
    >>> resolve(cfg, "reference", "fasta")
    PosixPath('/athena/.../data/ref/GRCh38_no_alt_analysis_set.fna')
    """
    node = cfg
    for k in keys:
        node = node[k]
    p = Path(node)
    if not p.is_absolute():
        p = cfg["root"] / p
    return p


def get_vcf_list(cfg: dict, platform: str) -> list[tuple[str, Path]]:
    """
    Return [(label, Path), ...] for a platform's VCF callsets.

    Parameters
    ----------
    platform : str
        One of 'novaseqx', 'ug100', 'ug100_ppmseq', 'ug100_ppmseq_duplex'.
    """
    key = f"{platform}_vcfs"
    entries = cfg.get(key, [])
    if not entries:
        return []
    return [(e["label"], resolve(cfg, key) if isinstance(entries, str)
             else cfg["root"] / e["vcf"])
            for e in entries]


# ============================================================
# Duplex read extraction (ppmSeq)
# ============================================================

def extract_duplex_reads(
    input_cram: Path,
    output_cram: Path,
    ref_fasta: Path,
    threads: int = 8,
) -> None:
    """
    Extract duplex reads from a ppmSeq CRAM.

    Duplex reads are identified by st and et tags both set to "MIXED".
    Per manuscript Methods: "we extracted duplex reads from the original
    CRAM files by selecting all reads where both st and et tag were set
    to 'MIXED'."

    Parameters
    ----------
    input_cram : Path to input ppmSeq CRAM
    output_cram : Path for output duplex-only CRAM
    ref_fasta : Path to reference FASTA
    threads : number of threads for samtools
    """
    output_cram.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "samtools", "view",
        "-@", str(threads),
        "-C",  # output CRAM
        "--reference", str(ref_fasta),
        "-e", '[st] == "MIXED" && [et] == "MIXED"',
        "-o", str(output_cram),
        str(input_cram),
    ]
    subprocess.run(cmd, check=True)

    # Index the output
    subprocess.run(["samtools", "index", str(output_cram)], check=True)


# ============================================================
# hap.py output parsing
# ============================================================

def read_pass_fp_fn(csv_path: Path) -> dict:
    """
    Parse hap.py summary.csv for PASS SNP and INDEL records.

    Returns
    -------
    dict : {'SNP': {'FP': float, 'FN': float, 'FP_FN_total': float},
            'INDEL': {...}}
    """
    out = {}
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    for vt in ("SNP", "INDEL"):
        match = [r for r in rows
                 if r.get("Type") == vt and r.get("Filter") == "PASS"]
        if not match:
            continue
        r = match[0]
        fp = float(r["QUERY.FP"])
        fn = float(r["TRUTH.FN"])
        out[vt] = {"FP": fp, "FN": fn, "FP_FN_total": fp + fn}
    return out


def read_extended_csv(ext_csv: Path) -> list[dict]:
    """Read hap.py extended.csv into a list of row dicts."""
    with ext_csv.open() as fh:
        return list(csv.DictReader(fh))


def extract_stratified_fpfn(
    rows: list[dict],
    subset: str,
    variant_type: str,
    filter_val: str = "PASS",
    subtype: str = "*",
) -> Optional[dict]:
    """
    Extract FP/FN/TP from hap.py extended CSV for a given Subset.

    Parameters
    ----------
    rows : list of dict from read_extended_csv()
    subset : str — e.g. '*', 'UG_HCR', 'UG_LCR', 'HP_10bp'
    variant_type : str — 'SNP' or 'INDEL'
    filter_val : str — 'PASS' or 'ALL'
    subtype : str — '*' by default

    Returns
    -------
    dict with keys FP, FN, TP, FP_FN_total, or None if not found.
    """
    match = [r for r in rows
             if r.get("Type") == variant_type
             and r.get("Filter") == filter_val
             and r.get("Subset") == subset
             and r.get("Subtype") == subtype]
    if not match:
        # Retry without Subtype constraint
        match = [r for r in rows
                 if r.get("Type") == variant_type
                 and r.get("Filter") == filter_val
                 and r.get("Subset") == subset]
    if not match:
        return None
    r = match[0]
    fp = float(r.get("QUERY.FP", r.get("FP", 0)))
    fn = float(r.get("TRUTH.FN", r.get("FN", 0)))
    tp = float(r.get("TRUTH.TP", r.get("TP", 0)))
    return {"FP": fp, "FN": fn, "TP": tp, "FP_FN_total": fp + fn}


# ============================================================
# Homopolymer stratum parsing (Figure 3)
# ============================================================

def parse_hp_length(subset_str: str) -> Optional[int]:
    """
    Extract numeric HP length from the Subset column.

    'HP_2bp'   -> 2
    'HP_gt20bp' -> 21  (plotted as '>20')

    Returns None if the string doesn't match.
    """
    if "gt20" in subset_str:
        return 21
    m = re.search(r"HP_(\d+)bp", subset_str)
    return int(m.group(1)) if m else None


def hp_label(hp_len: int) -> str:
    """Format HP length for display: 2-20 as-is, 21 as '>20'."""
    return str(hp_len) if hp_len <= 20 else ">20"


def load_hp_metrics(csv_path: Path, variant_type: str = "INDEL") -> dict:
    """
    Load hap.py extended.csv and extract per-HP-stratum metrics.

    Returns dict: hp_length -> {recall, precision, TP, FN, FP}
    """
    results = {}
    rows = read_extended_csv(csv_path)
    for row in rows:
        if row.get("Type") != variant_type:
            continue
        if row.get("Filter") != "PASS":
            continue
        subset = row.get("Subset", "")
        if not subset.startswith("HP_"):
            continue

        hp_len = parse_hp_length(subset)
        if hp_len is None:
            continue

        recall = row.get("METRIC.Recall", "")
        precision = row.get("METRIC.Precision", "")

        if recall and recall != "nan" and precision and precision != "nan":
            results[hp_len] = {
                "recall": float(recall),
                "precision": float(precision),
                "TP": float(row.get("TRUTH.TP", 0)),
                "FN": float(row.get("TRUTH.FN", 0)),
                "FP": float(row.get("QUERY.FP", 0)),
            }
    return results


# ============================================================
# samtools stats parsing (Figures 5 and 6B)
# ============================================================

def read_samtools_rl(rl_tsv: Path) -> tuple[list[int], list[int]]:
    """
    Parse samtools stats RL lines.

    Returns (read_lengths, counts).
    """
    lengths, counts = [], []
    with rl_tsv.open() as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0] == "RL":
                rl = int(parts[1])
                ct = int(parts[2])
                if ct > 0:
                    lengths.append(rl)
                    counts.append(ct)
    return lengths, counts


def compute_median_rl(lengths: list[int], counts: list[int]) -> float:
    """Compute median read length from a distribution."""
    total = sum(counts)
    cumsum = 0
    for rl, ct in zip(lengths, counts):
        cumsum += ct
        if cumsum >= total / 2:
            return float(rl)
    return 0.0


def read_samtools_gcd(gcd_tsv: Path) -> tuple[list[float], list[float]]:
    """
    Parse samtools stats GCD lines.

    Returns (gc_percentiles, median_normalized_depth).
    """
    gc_pcts, depths = [], []
    with gcd_tsv.open() as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 8 and parts[0] == "GCD":
                gc_pcts.append(float(parts[1]))
                depths.append(float(parts[5]))  # 50th percentile
    return gc_pcts, depths


def average_gcd(gcd_files: list[Path]) -> tuple[list[float], list[float]]:
    """Average normalized depth across replicates at each GC bin."""
    bin_depths = defaultdict(list)
    for f in gcd_files:
        gc, depth = read_samtools_gcd(f)
        for g, d in zip(gc, depth):
            bin_depths[g].append(d)
    sorted_bins = sorted(bin_depths.keys())
    avg_depth = [sum(bin_depths[g]) / len(bin_depths[g]) for g in sorted_bins]
    return sorted_bins, avg_depth


# ============================================================
# bedtools multiinter parsing (Figure 6A)
# ============================================================

def parse_multiinter(bed_path: Path) -> tuple[dict[int, int], int]:
    """
    Parse bedtools multiinter output.

    Returns
    -------
    counts_by_n : dict mapping n_shared -> total interval length
    n_samples : number of input samples
    """
    counts_by_n = {}
    n_samples = 0
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom\t"):
                parts = line.strip().split("\t")
                n_samples = len(parts) - 5
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            n_shared = int(parts[3])
            interval_len = int(parts[2]) - int(parts[1])
            counts_by_n[n_shared] = counts_by_n.get(n_shared, 0) + interval_len
    return counts_by_n, n_samples


def compute_reproducibility(counts_by_n: dict, n_samples: int) -> dict[int, float]:
    """
    Compute % reproducible variants for k = 2..n_samples.

    % at k = sum(counts where n >= k) / sum(all counts) * 100
    """
    total = sum(counts_by_n.values())
    if total == 0:
        return {}
    return {
        k: 100.0 * sum(v for n, v in counts_by_n.items() if n >= k) / total
        for k in range(2, n_samples + 1)
    }


# ============================================================
# Median across replicates
# ============================================================

def get_median_stats(label: str, per_rep_rows: list[dict]) -> dict:
    """
    Compute median FP and FN across replicates.

    Parameters
    ----------
    label : str — platform display name
    per_rep_rows : list of dicts from read_pass_fp_fn()

    Returns
    -------
    dict : {'SNP': {'platform': ..., 'FP': median, 'FN': median, ...}, 'INDEL': {...}}
    """
    out = {}
    for vt in ("SNP", "INDEL"):
        fps = [r[vt]["FP"] for r in per_rep_rows if vt in r]
        fns = [r[vt]["FN"] for r in per_rep_rows if vt in r]
        if not fps:
            continue
        med_fp = float(median(fps))
        med_fn = float(median(fns))
        out[vt] = {
            "platform": label,
            "variant_type": vt,
            "FP": med_fp,
            "FN": med_fn,
            "FP_FN_total": med_fp + med_fn,
        }
    return out


def compute_median_metrics(data_dict: dict) -> dict:
    """
    Compute median recall/precision across replicates at each stratum.

    Parameters
    ----------
    data_dict : dict mapping stratum_key -> list of {recall, precision, ...}

    Returns
    -------
    dict : stratum_key -> {recall, precision, recall_min, recall_max, ..., n_reps}
    """
    medians = {}
    for key in sorted(data_dict.keys()):
        reps = data_dict[key]
        medians[key] = {
            "recall": median([r["recall"] for r in reps]),
            "precision": median([r["precision"] for r in reps]),
            "recall_min": min(r["recall"] for r in reps),
            "recall_max": max(r["recall"] for r in reps),
            "precision_min": min(r["precision"] for r in reps),
            "precision_max": max(r["precision"] for r in reps),
            "n_reps": len(reps),
        }
    return medians


