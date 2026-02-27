#!/usr/bin/env python3
"""
scripts/figure_4/error_pileup.py

Empirical base-calling error rates by read position.

Per Supplementary Note 3:
  - MAPQ > 30, within NIST v4.2.1 high-confidence regions
  - Mask truth VCF variants
  - Traverse CIGAR: count mismatches, insertions, deletions
  - Stratify by read position and Q-score bin

Output: gzipped TSV per-position error counts.

Usage:
  python3 error_pileup.py \\
    --bam sample.bam --ref GRCh38.fna \\
    --regions NIST.bed --mask truth_positions.bed.gz \\
    --output sample_error_by_pos.tsv.gz \\
    --platform NovaSeqX_25B --label rep1
"""

import argparse
import gzip
import sys
from collections import defaultdict

import numpy as np
import pysam


def load_mask(mask_bed_gz):
    """Load truth variant mask as dict of chrom -> frozenset(positions)."""
    mask = defaultdict(set)
    with gzip.open(mask_bed_gz, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            for pos in range(start, end):
                mask[chrom].add(pos)
    return {c: frozenset(p) for c, p in mask.items()}


def load_regions(bed_path):
    """Load BED regions as list of (chrom, start, end)."""
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def process_bam(bam_path, ref_path, regions, mask, max_pos, platform):
    """Traverse reads and count errors by position in read."""
    max_q_threshold = 40 if ("NovaSeq" in platform or "NSX" in platform) else 35

    samfile = pysam.AlignmentFile(bam_path, reference_filename=ref_path)
    ref = pysam.FastaFile(ref_path)

    total = np.zeros(max_pos + 1, dtype=np.int64)
    mismatches = np.zeros(max_pos + 1, dtype=np.int64)
    insertions = np.zeros(max_pos + 1, dtype=np.int64)
    deletions = np.zeros(max_pos + 1, dtype=np.int64)
    total_maxQ = np.zeros(max_pos + 1, dtype=np.int64)
    mismatches_maxQ = np.zeros(max_pos + 1, dtype=np.int64)

    n_reads = 0
    for ri, (chrom, start, end) in enumerate(regions):
        if ri % 10000 == 0:
            print(f"  Region {ri}/{len(regions)} ({chrom}:{start}), reads: {n_reads:,}",
                  file=sys.stderr, flush=True)

        chrom_mask = mask.get(chrom, frozenset())

        for read in samfile.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality <= 30:
                continue

            n_reads += 1
            quals, seq = read.query_qualities, read.query_sequence
            if seq is None or quals is None:
                continue

            ref_start = read.reference_start
            try:
                ref_seq = ref.fetch(chrom, ref_start,
                                    ref_start + read.reference_length + 100)
            except (ValueError, KeyError):
                continue

            read_pos, ref_pos = 0, 0
            for op, length in read.cigartuples:
                if op == 0:  # M
                    for i in range(length):
                        rp = read_pos + i
                        if rp > max_pos:
                            break
                        genome_pos = ref_start + ref_pos + i
                        if genome_pos in chrom_mask or genome_pos < start or genome_pos >= end:
                            continue
                        total[rp] += 1
                        q = quals[rp]
                        if ref_pos + i < len(ref_seq):
                            if seq[rp].upper() != ref_seq[ref_pos + i].upper():
                                mismatches[rp] += 1
                                if q >= max_q_threshold:
                                    mismatches_maxQ[rp] += 1
                        if q >= max_q_threshold:
                            total_maxQ[rp] += 1
                    read_pos += length
                    ref_pos += length
                elif op == 1:  # I
                    rp = read_pos
                    if rp <= max_pos:
                        gp = ref_start + ref_pos
                        if gp not in chrom_mask and start <= gp < end:
                            insertions[rp] += length
                    read_pos += length
                elif op == 2:  # D
                    rp = read_pos
                    if rp <= max_pos:
                        gp = ref_start + ref_pos
                        if gp not in chrom_mask and start <= gp < end:
                            deletions[rp] += length
                    ref_pos += length
                elif op == 4:  # S
                    read_pos += length
                elif op == 3:  # N
                    ref_pos += length
                elif op in (7, 8):  # = or X
                    read_pos += length
                    ref_pos += length

    samfile.close()
    ref.close()
    print(f"  Processed {n_reads:,} reads", file=sys.stderr)

    return {
        "total": total, "mismatches": mismatches,
        "insertions": insertions, "deletions": deletions,
        "total_maxQ": total_maxQ, "mismatches_maxQ": mismatches_maxQ,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--max-position", type=int, default=350)
    args = parser.parse_args()

    print(f"Loading mask from {args.mask}...", file=sys.stderr)
    mask = load_mask(args.mask)

    print(f"Loading regions from {args.regions}...", file=sys.stderr)
    regions = load_regions(args.regions)

    print(f"Processing {args.bam}...", file=sys.stderr)
    counts = process_bam(args.bam, args.ref, regions, mask,
                         args.max_position, args.platform)

    with gzip.open(args.output, "wt") as fh:
        fh.write("platform\tlabel\tposition\ttotal_bases\tmismatches\t"
                 "insertions\tdeletions\ttotal_bases_maxQ\tmismatches_maxQ\n")
        for pos in range(args.max_position + 1):
            fh.write(f"{args.platform}\t{args.label}\t{pos}\t"
                     f"{counts['total'][pos]}\t{counts['mismatches'][pos]}\t"
                     f"{counts['insertions'][pos]}\t{counts['deletions'][pos]}\t"
                     f"{counts['total_maxQ'][pos]}\t{counts['mismatches_maxQ'][pos]}\n")

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

