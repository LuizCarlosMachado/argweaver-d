import os
import gzip
import argparse
import re

# Argument parsing
parser = argparse.ArgumentParser(description="Directly generate relabeled .newick files from .smc.gz")
parser.add_argument("positions_file", help="Path to positions file <chr> <pos>")
parser.add_argument("--min_iter", type=int, default=0, help="Minimum iteration number to consider (e.g., 1500)")
args = parser.parse_args()

# Base path
input_dir = "."
output_base_dir = "./target"

# Load labels
def get_labels_from_first_smc(target_dir="."):
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.endswith(".smc.gz"):
                with gzip.open(os.path.join(root, file), 'rt') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith("NAMES"):
                        return first_line.split()[1:]
    return []

labels = get_labels_from_first_smc()
if not labels:
    print("No labels found in any .smc.gz file inside ./target")
    exit(1)

# Build positions
position_labels = {}
with open(args.positions_file, "r") as f:
    for line in f:
        if line.strip():
            chrom, pos = line.strip().split()
            pos = int(pos)
            label = f"chr{chrom}_{pos}"
            position_labels[pos] = label

# Regex
nhx_pattern = re.compile(r'\[&&NHX:[^]]*\]')
node_pattern = re.compile(r'(\d+):')

def relabel_newick(newick):
    def replace_match(match):
        index = int(match.group(1))
        return labels[index] + ":" if 0 <= index < len(labels) else match.group(0)
    newick = re.sub(nhx_pattern, "", newick)
    return re.sub(node_pattern, replace_match, newick)

# Process
for filename in os.listdir(input_dir):
    if filename.endswith(".smc.gz"):
        match = re.search(r"\.(\d+)\.smc\.gz$", filename)
        if not match:
            continue
        iter_num = int(match.group(1))
        if iter_num < args.min_iter:
            continue

        input_path = os.path.join(input_dir, filename)
        with gzip.open(input_path, 'rt') as infile:
            current_regions = []
            for line in infile:
                if line.startswith("NAMES") or line.startswith("REGION"):
                    continue  # skip
                if line.startswith("TREE"):
                    parts = line.strip().split("\t")
                    tree_start = int(parts[1])
                    tree_end = int(parts[2])
                    for pos, label in position_labels.items():
                        if tree_start <= pos <= tree_end:
                            clean_newick = re.sub(r'^TREE\s+\d+\s+\d+\s+', '', line.strip())
                            relabeled = relabel_newick(clean_newick)

                            out_dir = os.path.join(output_base_dir, label)
                            os.makedirs(out_dir, exist_ok=True)
                            output_file = os.path.join(out_dir, filename.replace(".smc.gz", "_relabeled.newick"))

                            with open(output_file, 'w') as out:
                                out.write(relabeled + "\n")
                            print(f"[✓] Saved: {output_file}")

print(" All trees processed directly to relabeled.")
