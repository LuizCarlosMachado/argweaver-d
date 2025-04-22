#!/bin/bash

# Check for required argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 positions_file.txt"
    exit 1
fi

positions_file="$1"

# Read positions and build target directories
directories=()
while IFS=$'\t' read -r chr pos; do
    dir="target/chr${chr}_${pos}"
    directories+=("$dir")
done < "$positions_file"

# Process each existing directory
for dir in "${directories[@]}"; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        for file in "$dir"/*.smc.gz; do
            [ -e "$file" ] || continue  # Skip if no .smc.gz files
            output_file="${file/.smc.gz/.newick}"
            echo "Processing $file -> $output_file"
            zcat "$file" | sed 's/\[&&NHX:[^]]*\]//g' > "$output_file"
            echo "File saved: $output_file"
        done
    else
        echo "[!] Directory does not exist: $dir (skipped)"
    fi
done

echo "✅ All files have been processed!"

