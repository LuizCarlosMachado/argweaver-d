# -*- coding: utf-8 -*-
import gzip
import os
import argparse
import re

# Usage: python3 filter_smc.py chr${i}_positions.txt --min_iter 1500

# Command-line arguments
parser = argparse.ArgumentParser(description="Filter .smc.gz files by genomic positions and minimum iteration.")
parser.add_argument("positions_file", help="Path to the input file with positions (format: chr[TAB]pos)")
parser.add_argument("--min_iter", type=int, default=0, help="Minimum iteration number to consider (e.g., 1500)")
args = parser.parse_args()

output_base_dir = "./target"
input_dir = "."  # Directory with .smc.gz files

# Load positions from the provided file
position_labels = {}
with open(args.positions_file, "r") as f:
    for line in f:
        if line.strip() == "":
            continue
        chrom, pos = line.strip().split()
        pos = int(pos)
        label = f"chr{chrom}_{pos}"
        position_labels[pos] = label

# Process .smc.gz files based on minimum iteration
for filename in os.listdir(input_dir):
    if filename.endswith(".smc.gz"):
        match = re.search(r"\.(\d+)\.smc\.gz$", filename)
        if not match:
            continue
        iter_num = int(match.group(1))
        if iter_num < args.min_iter:
            continue  # skip files below the iteration threshold

        input_path = os.path.join(input_dir, filename)
        with gzip.open(input_path, 'rt') as infile:
            lines_by_position = {pos: [] for pos in position_labels}

            for line in infile:
                if line.startswith("NAMES") or line.startswith("REGION"):
                    for pos in lines_by_position:
                        lines_by_position[pos].append(line)
                elif line.startswith("TREE"):
                    parts = line.strip().split("\t")
                    tree_start = int(parts[1])
                    tree_end = int(parts[2])
                    for pos in lines_by_position:
                        if tree_start <= pos <= tree_end:
                            lines_by_position[pos].append(line)

        # Save filtered trees per position, only if data exists
        for pos, lines in lines_by_position.items():
            if len(lines) > 2:
                out_dir = os.path.join(output_base_dir, position_labels[pos])
                os.makedirs(out_dir, exist_ok=True)
                output_path = os.path.join(out_dir, filename)
                with gzip.open(output_path, 'wt') as outfile:
                    outfile.writelines(lines)
                print(f"[✓] Saved: {output_path} ({len(lines) - 2} trees)")
            else:
                print(f"[!] No trees for pos {pos} in {filename}, skipped.")

print("✅ Processing complete!")

