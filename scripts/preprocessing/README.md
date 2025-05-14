# Converting VCFs to `.sites` Files

The script (`vcf2sites.py`) converts phased and unphased VCFs into `.sites` format for use with ARGweaver and ARGweaver-D. It is designed for datasets that combine modern and ancient human genomes and assumes that all variants are **biallelic SNPs**.

---

## What the Script Does

* Parses genotypes from a VCF or VCF.gz file.
* Converts each sample’s genotype into a haplotype sequence.
* Assigns:

  * `REF/REF` → two reference alleles (e.g., "AA")
  * `ALT/ALT` → two alternate alleles (e.g., "GG")
  * 0/1, 1/0 Heterozygous → AG (e.g.)
  * Missing genotypes (`./.`) → `"NN"`
* Keeps only one line per position (removes duplicates).
* Filters out invariant sites (positions where all haplotypes are identical).
* Outputs a `.sites` file in the format required by ARGweaver.


### Note on Heterozygotes and Phasing

Although phased genotypes (such as 0|1) can indicate specific haplotypes, this script converts all heterozygous genotypes to "NN", regardless of phasing status. This is intentional. Since we work with mixed data — modern genomes (typically phased) and archaic genomes (typically unphased) — we chose a uniform approach to avoid introducing bias into the inference process. However, this behavior can be adjusted depending on the specific characteristics of your dataset.

---

## Usage

```bash
python3 vcf2sites.py input.vcf.gz output_prefix [chrom start end]
```

### Examples

**Full genome (no region filtering):**

```bash
python3 vcf2sites.py merged.vcf.gz chr22_sites
```

**Specific window (e.g., chr22:100000–200000):**

```bash
python3 vcf2sites.py merged.vcf.gz chr22_sites 22 100000 200000
```

---

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

## Notes

* Designed specifically for **biallelic SNPs only**. Multiallelic or indel variants should be filtered out beforehand (e.g., with `bcftools view -m2 -M2 -v snps`).
* The script ignores genotype phasing symbols (`/` vs `|`).
* Missing and heterozygous calls are encoded as `"NN"` to reflect uncertainty.

---




