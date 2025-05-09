# ARGweaver Tree Post-Processing Script

This script is designed for researchers working with **ARGweaverD** output who are interested in extracting and relabeling **local genealogies** at specific genomic positions.

---

## Purpose

This script enables users to:

* Extract genealogical trees corresponding to specific genomic positions of interest.
* Remove internal NHX metadata from Newick-formatted trees.
* Relabel node indices with sample names from the input VCF used during the ARGweaver run.

It is particularly useful for downstream analyses requiring correct sample labeling, such as allele frequency trajectory inference (e.g., CLUES2) or ARG-based selection scans.

---

## Input Requirements

1. **ARGweaverD output files** in `.smc.gz` format.

2. A plain-text file named `positions.txt` with tab-separated chromosome and position values, e.g.:

   ```
   2	35072772
   2	35385475
   2	35387871
   ```

3. (Optional) A minimum iteration threshold to exclude burn-in samples (default is 0).

---

## Usage

From within the directory containing your `.smc.gz` files, run:

```bash
python3 direct_relabel.py positions.txt --min_iter 1500
```

This will:

* Search for local trees covering each genomic position listed in `positions.txt`.
* Filter by the minimum iteration number, if specified.
* Save each cleaned and relabeled tree to the folder `target/chr<chrom>_<pos>/`.

---

## How Relabeling Works

ARGweaver internally represents leaf nodes with numeric indices. This script:

* Reads the original sample names from the first `NAMES` line found in the `.smc.gz` files.
* Replaces all numeric node labels in the extracted Newick trees with the corresponding sample names.

This ensures consistent and interpretable output across all trees and simplifies integration with downstream tools.

---

## Output

The script creates a subdirectory for each target position:

```
target/chr2_35072772/
├── <filename>_relabeled.newick
```

Each file contains a relabeled Newick tree corresponding to one iteration where the tree covered the specified site.

---




