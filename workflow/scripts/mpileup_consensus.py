#!/usr/bin/env python3
"""
mpileup_consensus.py
---------------------
Convert `samtools mpileup -aa` output to a FASTA consensus sequence.
Positions below --min-depth are written as 'N' (no reference bias).

Usage (inside Snakemake rule):
    samtools mpileup -aa -q 0 -Q 0 -f ref.fa aln.bam \\
      | python mpileup_consensus.py --min-depth 5 --sample SAMPLENAME \\
      > consensus.fa
"""

import sys
import re
import argparse
from collections import Counter


def parse_bases(bases_str: str, ref_base: str) -> list[str]:
    """Parse the pileup bases string and return a list of called bases."""
    bases = []
    i = 0
    ref_base = ref_base.upper()
    bases_str_upper = bases_str.upper()

    while i < len(bases_str):
        c = bases_str_upper[i]

        if c == "^":            # read-start marker followed by mapping-quality char
            i += 2
        elif c == "$":          # read-end marker
            i += 1
        elif c in "+-":         # insertion / deletion — skip the indel sequence
            i += 1
            num_str = ""
            while i < len(bases_str) and bases_str_upper[i].isdigit():
                num_str += bases_str_upper[i]
                i += 1
            i += int(num_str) if num_str else 0
        elif c in (".", ","):   # matches reference
            bases.append(ref_base)
            i += 1
        elif c in "ACGTN":      # explicit base call
            bases.append(c)
            i += 1
        elif c == "*":          # deletion placeholder
            bases.append("*")
            i += 1
        else:
            i += 1              # skip unknown characters

    return bases


def majority_base(bases: list[str]) -> str:
    """Return the most frequent non-deletion base, or 'N' if none."""
    real = [b for b in bases if b != "*"]
    if not real:
        return "N"
    return Counter(real).most_common(1)[0][0]


def main() -> None:
    parser = argparse.ArgumentParser(description="mpileup → FASTA consensus")
    parser.add_argument("--min-depth", type=int, default=5,
                        help="Minimum read depth; below this threshold output N")
    parser.add_argument("--sample", default=None,
                        help="Override FASTA header with this sample name")
    args = parser.parse_args()

    current_chrom: str | None = None
    sequence: list[str] = []

    def flush(chrom: str, seq: list[str]) -> None:
        header = args.sample if args.sample else chrom
        print(f">{header}")
        joined = "".join(seq)
        for i in range(0, len(joined), 60):
            print(joined[i:i + 60])

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue

        chrom    = parts[0]
        ref_base = parts[2]
        depth    = int(parts[3])
        raw_bases = parts[4]

        if chrom != current_chrom:
            if current_chrom is not None:
                flush(current_chrom, sequence)
            current_chrom = chrom
            sequence = []

        if depth < args.min_depth:
            sequence.append("N")
        else:
            bases = parse_bases(raw_bases, ref_base)
            sequence.append(majority_base(bases))

    if current_chrom is not None:
        flush(current_chrom, sequence)


if __name__ == "__main__":
    main()
