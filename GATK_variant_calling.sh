#!/bin/bash

# This is a script to call germline variants in a human WGS paired end read 2x100 bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

#======================= Tools =======================
# Using Ubuntu 22.04.2 through WSL2
#samtools version 1.16.1
#FastQC version 0.12.1
#multiqc version 1.14
#bwa version 0.7.17-r1188 (installed through conda)
#gatk4 version 4.4.0 (installed through conda)

#======================= Directories and Data =======================
# Prepare the directories
mkdir aligned_reads reads scripts results data supporting_files

# Download the data
cd reads
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
cd ..

#======================= Prep files =======================
cd supporting_files

# download reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# index ref -.fai file before running haplotype caller 
samtools faidx hg38.fa

# ref dict -.dict file before running haplotype caller
gatk CreateSequenceDictionary -R=hg38.fa -O=hg38.dict

# download known sites for BSQR (base quality score) from GATK resource bundle
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
cd ..

#======================= Variant Calling steps =======================
# Directories
ref="/mnt/c/Users/angel/Desktop/Rthings/GATK/supporting_files/hg38.fa"
known_sites="/mnt/c/Users/angel/Desktop/Rthings/GATK/supporting_files/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/mnt/c/Users/angel/Desktop/Rthings/GATK/aligned_reads"
reads="/mnt/c/Users/angel/Desktop/Rthings/GATK/reads" 
results="/mnt/c/Users/angel/Desktop/Rthings/GATK/results"
data="/mnt/c/Users/angel/Desktop/Rthings/GATK/data"

# -----------------------
# STEP 1: QC - Run FastQC
# -----------------------

fastqc -t 4 ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc -t 4 ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming, quality ok
# in this case, even if we had bad quality base calls we would not trim because 
# 1) the reads are already very short (100 bp) and further shortening would cause difficulties in mapping and
# 2) because bwa-mem performs softclipping

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

# BWA index reference
bwa index ${ref}

# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

samtools view SRR062634.paired.sam | less
samtools flagstat SRR062634.paired.sam

# ----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# ----------------------------------------

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634._sorted_dedup_reads.bam

samtools flagstat SRR062634._sorted_dedup_reads.bam

# ----------------------------------
# STEP 4: Base Quality Recalibration
# ----------------------------------

# build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634._sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

# apply the model to adjust base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634._sorted_dedup_reads.bam -R ${ref}  --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634._sorted_dedup_bqsr_reads.bam

# ------------------------------------------------
# STEP 5: Collect alignment and insert size metric
# ------------------------------------------------

gatk CollectAlignmentSummaryMetrics -R=${ref} -I=${aligned_reads}/SRR062634._sorted_dedup_bqsr_reads.bam -O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics -I=${aligned_reads}/SRR062634._sorted_dedup_bqsr_reads.bam -O=${aligned_reads}/insert_size_metrics.txt -H=${aligned_reads}/insert_size_histogram.pdf

# ---------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ---------------------------------------------

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634._sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf

# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
