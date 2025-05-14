# Converting VCFs to `.sites` Files

This repository provides two scripts for converting VCF or VCF.gz files into `.sites` format used by [ARGweaver](http://github.com/mdrasmus/argweaver) and [ARGweaver-D](https://github.com/CshlSiepelLab/argweaver-d). Both scripts assume **biallelic SNPs only** and are tailored for human genomic data (modern and ancient).

---

## Available Scripts

### ✅ `vcf2sites.py`
For **mixed datasets** (phased and unphased) — conservative handling.

- Best suited for datasets that include **ancient genomes** or low-quality data.
- Ignores phasing information to maintain uniform treatment across samples.

### ✅ `vcf2sites_phased.py`
For **fully phased** datasets.

- Assumes all genotypes are correctly phased (e.g., `0|1`, `1|0`).
- Splits genotypes into **haplotypes** using phasing information.
- Heterozygotes are properly resolved (e.g., `0|1` → `REF`, `ALT`).
- Best suited for **high-quality modern datasets**, such as from the 1000 Genomes Project.

---

## What Both Scripts Do

- Parse genotypes from `.vcf` or `.vcf.gz` files.
- Convert each diploid genotype into two haplotypes.
- Output a `.sites` file with:
  - A `NAMES` header line.
  - A `REGION` line with chromosome and coordinates.
  - One line per variant site: position and haplotype string.
- Remove invariant positions (where all haplotypes are identical).
- Skip positions with unresolvable or non-biallelic genotypes.
- Replace missing or malformed genotypes with `"NN"`.

---

## Input Requirements

Both scripts assume:

- **Biallelic SNPs only**  
  Before conversion, filter the input VCF with:

  ```bash
  bcftools view -m2 -M2 -v snps input.vcf.gz -Oz -o filtered.vcf.gz

## Output

A `.sites` file with the following structure:

```
NAMES   sample1_1 sample1_2 sample2_1 ...
REGION  22  100000  200000
101238  AAGG...
...
```

Each row represents a **biallelic SNP**, with haplotypes for each sample.

---




