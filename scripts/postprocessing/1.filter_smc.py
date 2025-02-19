# -*- coding: utf-8 -*-
import gzip
import os

# Load positions from positions.txt
positions_file = "positions.txt"
output_base_dir = "./target"

# Create the main target directory if it does not exist
os.makedirs(output_base_dir, exist_ok=True)

# Read positions from file
with open(positions_file, "r") as f:
    positions = {int(line.strip()): os.path.join(output_base_dir, line.strip()) for line in f}

# Create subdirectories for each position
for dir_path in positions.values():
    os.makedirs(dir_path, exist_ok=True)

# Process each .smc.gz file in the input directory
input_dir = "."  # Directory containing .smc.gz files
for filename in os.listdir(input_dir):
    if filename.endswith(".smc.gz"):
        input_path = os.path.join(input_dir, filename)

        with gzip.open(input_path, 'rt') as infile:
            keep_lines_by_position = {pos: [] for pos in positions}  # Store trees by position
            
            for line in infile:
                if line.startswith("NAMES") or line.startswith("REGION"):
                    for pos in positions:
                        keep_lines_by_position[pos].append(line)
                
                elif line.startswith("TREE"):
                    parts = line.strip().split("\t")
                    tree_start = int(parts[1])
                    tree_end = int(parts[2])

                    for pos, dir_path in positions.items():
                        if tree_start <= pos <= tree_end:
                            keep_lines_by_position[pos].append(line)

            # Save each file in the correct directory
            for pos, lines in keep_lines_by_position.items():
                if len(lines) > 2:
                    output_path = os.path.join(positions[pos], filename)
                    with gzip.open(output_path, 'wt') as outfile:
                        outfile.writelines(lines)
                    print(f"Filtered file saved: {output_path} with {len(lines) - 2} trees")
                else:
                    print(f"⚠️ No relevant trees found for position {pos} in {filename}, file ignored.")

print("Processing completed! Filtered trees are stored in './target/'.")

