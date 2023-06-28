#!/bin/bash

# This is a script to filter and annotate variants (continuation from GATK_variant_calling.sh)

#Directories 
ref="/mnt/c/Users/angel/Desktop/Rthings/Bioinformagician/GATK/supporting_files/hg38.fa"
results="/mnt/c/Users/angel/Desktop/Rthings/Bioinformagician/GATK/results"

# -----------------------
# Filter Variants - GATK4
# -----------------------

# Filter SNPs (following GATK4 best practices recommendations)
gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/raw_snps.vcf \
    -O ${results}/filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter" \

# Filter INDELS
gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/raw_indels.vcf \
    -O ${results}/filtered_indels.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter" \

# Select Variants that PASS the filters (those tagged "filter-name" above)
gatk SelectVariants --exclude-filtered -V ${results}/filtered_snps.vcf -O ${results}/analysis-ready-snps.vcf
gatk SelectVariants --exclude-filtered -V ${results}/filtered_indels.vcf -O ${results}/analysis-ready-indels.vcf 

# Select Variants that PASS the filters (those tagged "genotype-filter" above)
cd results
# first check that the grep command got the desired result
cat analysis-ready-snps.vcf | grep -v -E "DP_filter|GQ_filter" | less
# now save grep output to another file
cat analysis-ready-snps.vcf | grep -v -E "DP_filter|GQ_filter" > analysis-ready-snps-filteredGT.vcf
cat analysis-ready-indels.vcf | grep -v -E "DP_filter|GQ_filter" > analysis-ready-indels-filteredGT.vcf

# ------------------------------------
# Annotate Variants - GATK4 Funcotator
# ------------------------------------

# Download funcotator's germline data sources (if first time using it)
# ./GATK FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download

# Annotate SNPs using Funcotator
gatk Funcotator \
    --variant ${results}/analysis-ready-snps-filteredGT.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path path/to/data/sources \
    --output ${results}/analysis-ready-snps-filteredGT-funcotated.vcf \
    --output-file-format VCF 

# Annotate INDELs using Funcotator
gatk Funcotator \
    --variant ${results}/analysis-ready-indels-filteredGT.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path  path/to/data/sources \
    --output ${results}/analysis-ready-indels-filteredGT-funcotated.vcf \
    --output-file-format VCF 

# Extract fields from a VCF file to a tab-delimited file
gatk VariantsToTable \
    -V ${results}/analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ${results}/output_snps.table

cd results
cat output_snps.table | less
cat analysis-ready-snps-filteredGT-funcotated.vcf | less
# select the funcotator field only
cat analysis-ready-snps-filteredGT-funcotated.vcf | grep " Funcotation fields are: " | less
# replace | (pipe operator) with tab
cat analysis-ready-snps-filteredGT-funcotated.vcf | grep " Funcotation fields are: " | sed 's/|/\t/g' | less
cat analysis-ready-snps-filteredGT-funcotated.vcf | grep " Funcotation fields are: " | sed 's/|/\t/g' > output_curated_variants.txt
# extract the values in funcotation column (5th field), grep gene of interest e.g. "NBPF1" 
# (so we pick only the rows from the column funcotation that have to do with this gene)
# and then replace pipe operator with tab
cat output_snps.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' | less
# append these values (>>) to the output_curated_variants.txt file
cat output_snps.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' >> output_curated_variants.txt

# In the end we get a txt file in tabular form with information about the SNPs of a specific gene (e.g. if they are in introns,
# exons, 5 prime flank etc) which is very helpful.
