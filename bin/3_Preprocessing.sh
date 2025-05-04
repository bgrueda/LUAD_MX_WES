#!/bin/bash

# Samtools (v.1.10) have some tools to manage the sam and bam files and do the depth and coverage analysis.
# GATK (v.4.1.9.0) and Picard (v.2.22.8) to mark duplicates and do the base recalibration.

Res="./../results"
Ref="/labs/lgfc/hg19/"

for i in $Res*sort.bam
do
	time gatk MarkDuplicates \
        I=$i O=${i%_sort.bam}'_md.bam' \
        M=${i%_sort.bam}_metrics.txt 
	time gatk BaseRecalibrator \
        -I ${i%_sort.bam}'_md.bam' \
		-R $Ref'Homo_sapiens_assembly19.fasta'  \
		--known-sites $Ref'Homo_sapiens_assembly19.dbsnp138.vcf' \
		--known-sites $Ref'Homo_sapiens_assembly19.known_indels.vcf' \
		--known-sites $Ref'Mills_and_1000G_gold_standard.indels.b37.sites.vcf' \
		-O ${i%_sort.bam}.table
	time gatk ApplyBQSR -R $Ref'Homo_sapiens_assembly19.fasta' \
        -I ${i%_sort.bam}'_md.bam' \
		--bqsr-recal-file ${i%_sort.bam}.table \
		-O ${i%_sort.bam}_br.bam
	samtools index -b ${i%_sort.bam}'_br.bam'
	
done

for i in $Res*br.bam
do
    samtools view -b -q 20 $i > ${i%_br.bam}'_ready.bam'
    samtools index -b ${i%_br.bam}'_ready.bam'
done