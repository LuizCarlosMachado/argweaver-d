#!/bin/bash

# Read positions from positions.txt
directories=()
while IFS= read -r pos; do
    directories+=("target/$pos")
done < positions.txt

# Process each directory
for dir in "${directories[@]}"; do
    echo "📂 Processing directory: $dir"

    for file in "$dir"/*.smc.gz; do
        output_file="${file/.smc.gz/.newick}"

        echo "🛠️  Processing $file -> $output_file"
        zcat "$file" | sed 's/\[&&NHX:[^]]*\]//g' > "$output_file"

        echo "File saved: $output_file"
    done
done

echo "All files have been processed!"

