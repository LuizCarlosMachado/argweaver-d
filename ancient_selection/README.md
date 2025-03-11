```markdown
# Analysis of Highly Frequent Variants in Modern Humans and Their Absence in Ancient DNA

## **Overview**
This analysis aims to identify genomic regions that are **highly frequent in modern humans (`AF ≥ 0.95`)** but **absent in ancient DNA (aDNA)**. Such regions are potential candidates for **ancient positive selection** that contributed to the success of our species.

## **Workflow**

### **1. Filtering High-Frequency Variants in Modern Humans (1000 Genomes Project)**
- Extracted **biallelic SNPs** with **allele frequency (`AF`) ≥ 0.95** from modern human VCFs:
  ```bash
  bcftools view -i 'INFO/AF>=0.95' -m2 -M2 -v snps ALL.chr*.vcf.gz -Oz -o highFrequency/chr*_biallelic.vcf.gz
  ```
- Saved the positions of these variants for later filtering:
  ```bash
  cut -d',' -f1,2 --output-delimiter=$'\t' biallelic_high_freq_pos.csv > regions.txt
  ```

### **2. Subsetting Ancient VCFs (Altai, Chagyrskaya, Denisova, Vindija)**
- Filtered ancient DNA VCFs to retain **only the high-frequency modern human variants**:
  ```bash
  bcftools view -R regions.txt -Oz -o Altai_subset.vcf.gz Altai.vcf.gz
  bcftools view -R regions.txt -Oz -o Chagyrskaya_subset.vcf.gz Chagyrskaya.vcf.gz
  bcftools view -R regions.txt -Oz -o Denisova_subset.vcf.gz Denisova.vcf.gz
  bcftools view -R regions.txt -Oz -o Vindija_subset.vcf.gz Vindija.vcf.gz
  ```

### **3. Merging Ancient VCFs**
- Combined all ancient VCFs into a single file:
  ```bash
  bcftools merge -m all -Oz -o ancient_merged.vcf.gz Altai_subset.vcf.gz Chagyrskaya_subset.vcf.gz Denisova_subset.vcf.gz Vindija_subset.vcf.gz
  bcftools index ancient_merged.vcf.gz
  ```

### **4. Extracting Genotypes**
- Removed all additional fields, keeping only genotype (`GT`) data:
  ```bash
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' ancient_merged.vcf.gz > ancient_merged_genotypes.tsv
  ```

### **5. Analyzing Genotype Combinations**
- Counted occurrences of unique genotype combinations in ancient samples:
  ```python
  import pandas as pd

  df = pd.read_csv("ancient_merged_genotypes.tsv", sep="\t", header=None)
  df = df.iloc[:, 1:]  # Remove the position column
  df["Genotype Combination"] = df.apply(lambda row: "|".join(row.astype(str)), axis=1)
  genotype_counts = df["Genotype Combination"].value_counts().reset_index()
  genotype_counts.columns = ["Genotype Combination", "Count"]
  genotype_counts.to_csv("genotype_combination_counts.csv", index=False)
  print(genotype_counts.head(10))
  ```
