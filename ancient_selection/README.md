# **Identifying Divergent Variants Between Modern Humans and Archaic Genomes**

## **Objective**
The goal of this pipeline is to identify genomic positions where modern humans and archaic hominins (Neanderthals and Denisovans) differ significantly. We focus on variants that are **highly frequent in modern human populations (1000 Genomes Project)** but are **absent or fixed differently in archaic genomes**.

To ensure we capture all relevant variants, we apply two complementary approaches:
1. **Variants fixed in modern humans (1/1) but absent (0/0) in archaic genomes.**
2. **Variants fixed in archaic genomes (1/1) but absent (0/0) in modern humans.**

These two categories may reveal functionally significant mutations, including **adaptive changes, reference genome errors, or population-specific selection events**.

---

# **Part 1: Identifying Variants Fixed as 1/1 in Modern Humans and Absent (0/0) in Archaic Genomes**
### **Rationale**
To detect variants that became highly frequent in modern humans but were absent in archaic populations, we:
- Identify positions where **≥95% of modern humans are 1/1**.
- Check if these positions exist in archaic genomes.
- Filter for cases where archaic genomes are **fixed as 0/0**.
- Ensure that both modern and archaic datasets contain the **exact same positions** before merging.

### **Step-by-Step Execution**
#### **1. Filtering High-Frequency ALT Variants in Modern Humans**
```bash
bcftools view -i 'INFO/AF >= 0.95' -m2 -M2 -v snps ALL.chr*.vcf.gz -Oz -o high_AF_modern_chr*.vcf.gz
```

#### **2. Filtering Low-Frequency REF Variants in Modern Humans**
```bash
bcftools view -i 'INFO/AF <= 0.05 & INFO/AF > 0' -m2 -M2 -v snps ALL.chr*.vcf.gz -Oz -o low_AF_modern_chr*.vcf.gz
```

Both `high_AF_modern_chr` and `low_AF_modern_chr` were processed identically.

---

#### **3. Index and Concatenate High AF VCFs**
```bash
for file in high_AF_modern_chr*.vcf.gz; do
    bcftools index -t "$file"
done

bcftools concat -Oz -o concatenated_high_AF_modern.vcf.gz $(ls high_AF_modern_chr*.vcf.gz | sort -V)
```

---

#### **4. Extract Positions for High-Frequency Variants in Modern Humans**
```bash
bcftools query -f '%CHROM\t%POS\n' concatenated_high_AF_modern.vcf.gz > high_AF_positions.txt
```

#### **5. Subset Archaic Genomes for These Positions**
```bash
bcftools view -R high_AF_positions.txt -Oz -o Altai_high_AF.vcf.gz Altai.vcf.gz
bcftools view -R high_AF_positions.txt -Oz -o Chagyrskaya_high_AF.vcf.gz Chagyrskaya.vcf.gz
bcftools view -R high_AF_positions.txt -Oz -o Denisova_high_AF.vcf.gz Denisova.vcf.gz
bcftools view -R high_AF_positions.txt -Oz -o Vindija_high_AF.vcf.gz Vindija.vcf.gz
```

---

#### **6. Merge Archaic VCFs**
```bash
bcftools merge -m all -Oz -o ancient_high_AF_merged.vcf.gz Altai_high_AF.vcf.gz Chagyrskaya_high_AF.vcf.gz Denisova_high_AF.vcf.gz Vindija_high_AF.vcf.gz
bcftools index ancient_high_AF_merged.vcf.gz
```

---

#### **7. Extract Positions Fixed as 0/0 in Archaic Genomes**
```bash
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ancient_high_AF_merged.vcf.gz > ancient_high_AF_merged.tsv
```
```bash
awk '$3=="0/0" && $4=="0/0" && $5=="0/0" && $6=="0/0" {print $1"\t"$2}' ancient_high_AF_merged.tsv > ancient_absent_ref_positions.txt
```

---

#### **8. Subset Modern and Archaic VCFs**
```bash
bcftools view -R ancient_absent_ref_positions.txt -Oz -o subset_ancient_high_AF.vcf.gz ancient_high_AF_merged.vcf.gz
bcftools view -R ancient_absent_ref_positions.txt -Oz -o subset_modern_high_AF.vcf.gz concatenated_high_AF_modern.vcf.gz
```

---

#### **9. Merge Final Filtered VCFs**
```bash
bcftools merge -m all -Oz -o final_merged_1_1_modern_0_0_ancient.vcf.gz subset_ancient_high_AF.vcf.gz subset_modern_high_AF.vcf.gz
```

---

# **Part 2: Identifying Variants Fixed as 1/1 in Archaic Genomes and Absent (0/0) in Modern Humans**
### **Rationale**
Instead of starting with modern humans, this pipeline **first cleans the archaic genomes** to ensure we work only with variants that are **1/1 across all archaic individuals**. This avoids unnecessary processing of **0/0 and missing (`./.`) genotypes**.

---

### **Step-by-Step Execution**
#### **1. Filter Archaic VCFs to Keep Only 1/1 Variants**
```bash
bcftools view -i 'GT="1/1"' -Oz -o Altai_filtered.vcf.gz Altai.vcf.gz
bcftools view -i 'GT="1/1"' -Oz -o Chagyrskaya_filtered.vcf.gz Chagyrskaya.vcf.gz
bcftools view -i 'GT="1/1"' -Oz -o Denisova_filtered.vcf.gz Denisova.vcf.gz
bcftools view -i 'GT="1/1"' -Oz -o Vindija_filtered.vcf.gz Vindija.vcf.gz
```

#### **2. Merge Filtered Archaic VCFs**
```bash
bcftools merge -m all -Oz -o ancient_filtered_merged.vcf.gz Altai_filtered.vcf.gz Chagyrskaya_filtered.vcf.gz Denisova_filtered.vcf.gz Vindija_filtered.vcf.gz
```

#### **3. Remove Positions with Missing Genotypes (`./.`)**
```bash
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ancient_filtered_merged.vcf.gz > ancient_filtered_merged.tsv
awk '$3=="./." || $4=="./." || $5=="./." || $6=="./." {print $1"\t"$2}' ancient_filtered_merged.tsv > positions_with_missing.txt
bcftools view -T ^positions_with_missing.txt -Oz -o ancient_final_filtered.vcf.gz ancient_filtered_merged.vcf.gz
```

#### **4. Extract Positions from Modern VCF**
```bash
bcftools query -f '%CHROM\t%POS\n' concatenated_low_AF_modern.vcf.gz > low_AF_positions.txt
```

#### **5. Subset Archaic VCF for These Positions**
```bash
bcftools view -R low_AF_positions.txt -Oz -o subset_ancient_low_AF.vcf.gz ancient_final_filtered.vcf.gz
```

#### **6. Extract Positions from the New Archaic Subset**
```bash
bcftools query -f '%CHROM\t%POS\n' subset_ancient_low_AF.vcf.gz > final_matched_positions.txt
```

#### **7. Subset Modern VCF for These Positions**
```bash
bcftools view -R final_matched_positions.txt -Oz -o subset_modern_low_AF.vcf.gz concatenated_low_AF_modern.vcf.gz
```

#### **8. Merge Final Filtered VCFs**
```bash
bcftools merge -m all -Oz -o final_merged_0_0_modern_1_1_ancient.vcf.gz subset_ancient_low_AF.vcf.gz subset_modern_low_AF.vcf.gz
```

---

## **Conclusion**
This two-part pipeline efficiently identifies variants **fixed in modern or archaic humans**, reducing noise and ensuring a clean dataset for evolutionary and functional analysis.
