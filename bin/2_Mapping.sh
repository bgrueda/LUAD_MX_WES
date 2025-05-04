#!/bin/bash

# BWA (v.0.7.12-r1039) align the reads to a reference genome
# You need to download the reference genome version that you prefer
# In this case,we are going to use: hg19
# It can be downloaded from: https://data.broadinstitute.org/snowman/hg19/

# Short routes to the directories
Data="/labs/lgfc/Fasta/"
Res="./../results"s
Ref="/labs/lgfc/hg19/" # Directory with reference genome files 

for i in $Data*1.fq.gz
do
	name="${i##*/}"
	nombre="${name%%'_'*}"
	bwa mem -R '@RG\tID:'$nombre'\tSM:'$nombre'\tPL:Illumina' \
        $Ref'Homo_sapiens_assembly19.fasta' \
        $i ${i%1.fq.gz}2.fq.gz > $Res$nombre'.sam'  
	samtools view -b -o $Res$nombre'.bam'  $Res$nombre'.sam'
	samtools sort $Res$nombre'.bam' -o $Res$nombre'_sort.bam'
	
done