## **Identifying High-Frequency Variants in Modern Humans and Their Absence in Archaic Genomes**

### **Objective**
The goal of this pipeline is to identify genomic positions in modern human populations from the **1000 Genomes Project** where variants are at very **high frequency (AF ≥ 95%)** but are **absent in archaic genomes** (Neanderthals and Denisovans).

To ensure we capture all relevant variants, we split the filtering into two parts:
- **High-frequency alternative alleles (AF ≥ 95%)** in modern humans.
- **Low-frequency reference alleles (AF ≤ 5%)** in modern humans.  
  This second step is important to detect potential reference genome errors where an allele may be 0/0 in modern humans and 1/1 in archaic genomes. These sites could either be fixed in archaic humans or result from a misrepresentation in the reference genome, making them worth investigating for selection.

---

## **Pipeline Overview**

### **High AF in Modern Humans**
- Identify **positions where ≥95% of modern humans are 1/1**.
- Check which of these variants are present in archaic genomes.
- Identify which of these variants are **0/0 in archaic genomes**.
- Perform subset filtering to retain only positions **1/1 in modern humans and 0/0 in archaic humans**.

### **Low AF in Modern Humans**
- Identify **positions where ≤5% of modern humans carry an alternative allele (mostly 0/0)**.
- Check which of these variants are present in archaic genomes.
- Identify which of these variants are **1/1 in archaic genomes**.
- Perform subset filtering to retain only positions **0/0 in modern humans and 1/1 in archaic humans**.

---

## **Step-by-Step Execution**

### **1. Filtering High-Frequency ALT Variants in Modern Humans**
```bash
bcftools view -i 'INFO/AF >= 0.95' -m2 -M2 -v snps ALL.chr*.vcf.gz -Oz -o high_AF_modern_chr*.vcf.gz
```

### **2. Filtering Low-Frequency REF Variants in Modern Humans**
```bash
bcftools view -i 'INFO/AF <= 0.05 & INFO/AF > 0' -m2 -M2 -v snps ALL.chr*.vcf.gz -Oz -o low_AF_modern_chr*.vcf.gz
```

For consistency, all steps applied to `high_AF_modern_chr` were also applied to `low_AF_modern_chr`.

---

## **3. Index and Concatenate High AF VCFs**
```bash
for file in low_AF_modern_chr*.vcf.gz; do
    bcftools index -t "$file"
done

bcftools concat -Oz -o concatenated_high_AF_modern.vcf.gz $(ls high_AF_modern_chr*.vcf.gz | sort -V)
```

---

## **4. Extract Positions for High-Frequency Variants in Modern Humans**
```bash
bcftools query -f '%CHROM\t%POS\n' concatenated_high_AF_modern.vcf.gz > high_AF_positions.txt
```

Using these positions, we subset the archaic genomes:
```bash
bcftools view -R high_AF_positions.txt -Oz -o Altai_high_AF.vcf.gz Altai.vcf.gz
bcftools view -R high_AF_positions.txt -Oz -o Chagyrskaya_high_AF.vcf.gz Chagyrskaya.vcf.gz
bcftools view -R high_AF_positions.txt -Oz -o Denisova_high_AF.vcf.gz Denisova.vcf.gz
bcftools view -R high_AF_positions.txt -Oz -o Vindija_high_AF.vcf.gz Vindija.vcf.gz
```

---

## **5. Merge Archaic VCFs**
```bash
bcftools merge -m all -Oz -o ancient_high_AF_merged.vcf.gz Altai_high_AF.vcf.gz Chagyrskaya_high_AF.vcf.gz Denisova_high_AF.vcf.gz Vindija_high_AF.vcf.gz
bcftools index ancient_high_AF_merged.vcf.gz
```

At this point, `ancient_high_AF_merged.vcf.gz` contains genomic positions that are **highly frequent in modern populations**, but we need to refine the dataset to retain only positions where modern humans are **1/1** and archaic genomes are **0/0**.

---

## **6. Extract Absent Positions as 0/0 in Archaic Genomes**
```bash
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ancient_high_AF_merged.vcf.gz > ancient_high_AF_merged.tsv
```
Using **awk**, we filter for positions where all four archaic genomes are **0/0**:
```bash
awk '$3=="0/0" && $4=="0/0" && $5=="0/0" && $6=="0/0" {print $1"\t"$2}' ancient_high_AF_merged.tsv > ancient_absent_ref_positions.txt
```

---

## **7. Subset Modern and Archaic VCFs**
Filter the archaic VCF to retain only 0/0 positions:
```bash
bcftools view -R ancient_absent_ref_positions.txt -Oz -o subset_ancient_high_AF.vcf.gz ancient_high_AF_merged.vcf.gz
```
Filter the modern VCF for the same positions:
```bash
bcftools view -R ancient_absent_ref_positions.txt -Oz -o subset_modern_high_AF.vcf.gz concatenated_high_AF_modern.vcf.gz
```

---

## **8. Merge Final Filtered VCFs**
```bash
bcftools merge -m all -Oz -o final_merged_1_1_modern_0_0_ancient.vcf.gz subset_ancient_high_AF.vcf.gz subset_modern_high_AF.vcf.gz
```
This final merged VCF (`final_merged_1_1_modern_0_0_ancient.vcf.gz`) contains only **positions where modern humans are fixed for the alternate allele (1/1) and archaic genomes are absent for the reference allele (0/0)**.

---

## **Conclusion**
This pipeline successfully identifies genomic positions where:
- **Modern humans (1000 Genomes Project) are fixed for the alternate allele (1/1)**.
- **Archaic genomes (Neanderthals and Denisovans) are fixed for the reference allele (0/0)**.
