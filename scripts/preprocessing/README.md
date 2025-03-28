# Generating ARGweaver-compatible .sites files from merged VCFs

This document describes the rationale, preprocessing steps, and the unified Python script used to convert merged VCF files (containing modern and ancient genomes) into `.sites` format compatible with ARGweaver. The approach ensures proper handling of genotypes, phasing, and invariant site filtering.

---

## Context

When preparing data for ancestral recombination graph (ARG) inference using ARGweaver, it is often necessary to merge phased VCFs from modern populations (e.g., 1000 Genomes Project) with pseudo-haploid or unphased VCFs from ancient genomes (e.g., Neanderthals, Denisovans). These merged datasets present unique challenges:

- Ancient VCFs may lack standard `ALT` fields.
- Modern VCFs contain only variant sites, while ancient VCFs often represent dense genomic sequences.
- Ancient genotypes are typically unphased, while modern ones are phased.

To standardize and densify the merged VCF, we generate a `.sites` file with:

- All informative positions represented;
- Correct haplotype encoding per individual;
- Optional exclusion of non-human outgroups (e.g., panTro4);
- Filtering of invariant sites to reduce file size and focus inference on polymorphic regions.

---

## Genotype Handling Logic

Each genotype (GT field) is parsed and translated into a pair of haplotype bases using the following rules:

| Genotype   | Meaning                 | Haplotype Output |          |
| ---------- | ----------------------- | ---------------- | -------- |
| 0/0 or 0   | 0                       | Homozygous REF   | REF, REF |
| 1/1 or 1   | 1                       | Homozygous ALT   | ALT, ALT |
| 0/1 or 1/0 | Heterozygous (unphased) | N, N             |          |
| ./., .     | Missing                 | N, N             |          |

- All phasing symbols (`|`) are treated the same as `/` for consistency.
- Only the first two fields of each genotype are used.
- Non-numeric or out-of-range allele indices are mapped to `N`.
- ALT alleles are extracted only if present; otherwise, only REF is used.

---

## Script Behavior

The Python script performs the following steps:

1. Reads the input `.vcf` or `.vcf.gz` file within a given chromosome window.
2. Converts each genotype into two haplotypes.
3. Constructs `.sites` lines with the position and the concatenated haplotypes.
4. Removes invariant sites (where all haplotypes have the same base).
5. Writes the output to a `.sites` file.

The script also supports an optional `--exclude-last` flag to remove the last individual (commonly used for outgroups like panTro4) from the analysis.

---

## Usage

```bash
python3 vcf2filtered_sites.py input.vcf.gz output_prefix chrom start end [--exclude-last]
```

### Example

```bash
python3 vcf2filtered_sites.py merged_chr22.vcf.gz chr22_sites 22 19000000 20000000 --exclude-last
```

---

## Output

The final output is a valid `.sites` file with the following structure:

```text
NAMES   sample1_1   sample1_2   ...
REGION  22  19000000    20000000
19173496    CCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
...
```

Only variant positions are retained after internal filtering.

---

## Script Location

The full script is available in the `script/` directory under the name:

```
vcf2filtered_sites.py
```

It is self-contained, and no external libraries beyond Python's standard library are required.

---

## Notes

- This pipeline assumes that merged VCFs have been harmonized in terms of REF/ALT alleles and contain no structural variants.
- For ancient genomes with missing genotypes (`./.`), the script assigns `N` to both haplotypes.
- Invariant positions are filtered **after** translation into haplotypes, to ensure accuracy even in the presence of missing data.

For questions or improvements, please submit an issue or pull request.

