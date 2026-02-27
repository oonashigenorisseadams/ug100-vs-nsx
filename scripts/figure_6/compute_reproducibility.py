#!/usr/bin/env python3
"""
scripts/figure_6/compute_reproducibility.py

Reads bedtools multiinter output and computes % reproducible variants
for k=2..n replicates. Called by run_reproducibility.sh.
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.benchmark_utils import parse_multiinter, compute_reproducibility


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiinter-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mi_dir = Path(args.multiinter_dir)
    rows = []

    for platform in ["NSX", "UG"]:
        for vtype in ["snps", "indels"]:
            for suffix, region in [
                ("_multiinter.bed", "Whole Genome"),
                ("_multiinter_NIST.bed", "NIST v4.2.1"),
            ]:
                bed = mi_dir / f"{platform}_{vtype}{suffix}"
                if not bed.exists():
                    print(f"  WARN: Missing {bed}")
                    continue

                counts, n_samples = parse_multiinter(bed)
                if n_samples == 0:
                    continue

                repro = compute_reproducibility(counts, n_samples)
                for k, pct in repro.items():
                    vt_label = "SNP" if "snp" in vtype else "INDEL"
                    rows.append({
                        "platform": "NovaSeqX 25B" if platform == "NSX" else "UG 100",
                        "variant_type": vt_label,
                        "region": region,
                        "n_reps_required": k,
                        "n_total_reps": n_samples,
                        "pct_reproducible": f"{pct:.2f}",
                    })
                    print(f"  {platform} {vtype} {region}: k={k} -> {pct:.1f}%")

    out_path = Path(args.output)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "platform", "variant_type", "region",
            "n_reps_required", "n_total_reps", "pct_reproducible",
        ])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
