# -*- coding: utf-8 -*-
import os
import re
import gzip

# Read labels directly from the first .smc.gz file in ./target
def get_labels_from_smc():
    target_dir = "./target"
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.endswith(".smc.gz"):
                with gzip.open(os.path.join(root, file), 'rt') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith("NAMES"):
                        return first_line.split()[1:]  # Ignore "NAMES" and return labels
    return []

# Load labels dynamically
labels = get_labels_from_smc()
if not labels:
    print("No labels found in .smc.gz files. Check your files!")
    exit()

# Read positions from positions.txt
with open("positions.txt", "r") as f:
    directories = [f"target/{line.strip()}" for line in f]

# Regular expressions for cleaning NHX metadata and replacing node numbers
nhx_pattern = re.compile(r'\[&&NHX:[^]]*\]')
node_pattern = re.compile(r'(\d+):')

# Function to replace indices with correct labels
def relabel_newick(newick):
    def replace_match(match):
        index = int(match.group(1))
        return labels[index] + ":" if 0 <= index < len(labels) else match.group(0)

    newick = re.sub(nhx_pattern, "", newick)
    return re.sub(node_pattern, replace_match, newick)

# Process each directory
for dir_path in directories:
    print(f"Processing directory: {dir_path}")

    for file_name in os.listdir(dir_path):
        if file_name.endswith(".newick"):
            file_path = os.path.join(dir_path, file_name)

            print(f"Processing {file_path}")

            with open(file_path, "r") as file:
                lines = file.readlines()

            tree_lines = [line for line in lines if line.startswith("TREE")]
            if not tree_lines:
                print(f"⚠️  No tree found in {file_path}")
                continue

            clean_newick = re.sub(r'^TREE\s+\d+\s+\d+\s+', '', tree_lines[0]).strip()
            relabeled_newick = relabel_newick(clean_newick)

            output_file_path = os.path.join(dir_path, file_name.replace(".newick", "_relabeled.newick"))
            with open(output_file_path, "w") as output_file:
                output_file.write(relabeled_newick)

            print(f" File saved: {output_file_path}")

print("All trees have been processed and corrected!")

