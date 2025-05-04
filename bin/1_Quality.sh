#!/bin/bash

# fastQC is for the quality control of raw data.

# First, be sure about have all the software from the previous requirements installed.
# Now, put all the fastq files with the raw data in the same directory.

# Create a new one to put in the results.
mkdir ./../results/1_Quality

# Run FastQC for the analysis
for i in ~/../data/*.fastq;
do fastqc $i;
done

# Move your data to a different directory
mv *.zip *.html ./../results/1_Quality/

# You need to examinate every one individually as is described on:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

#--------------------------------------------------------------------------------------------------------------------------
# Trimmomatic is the tool that allows you to do a quality filtering of the raw data.

# This part coul change depending of the results of the fastQC
# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 -> Remove adapters
# LEADING:3 -> Remove leading low quality or N bases (below quality 3)
# TRAILING:3 -> Remove trailing low quality or N bases (below quality 3)
# SLIDINGWINDOW:4:15 -> Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
# MINLEN:50 -> Drop reads below the 50 bases long.
# CROP: 150 -> The number of bases to keep, from the start of the read.
# HEADCROP:15 -> The number of bases to remove from the start of the read.
#--------------------------------------------------------------------------------------------------------------------------

# OPTINAL. If the raw data quality is insufficient, trimming steps can be applied to enhance it and achieve the required quality. 

for i in ./../data/*R1.fastq;
do java -jar /route/trimmomatic-0.39.jar PE -phred33 $i ${i%?.fastq}2.fastq \
    trimmed_$i \
    unpaired_$1 \
    trimmed${i%?.fastq}2.fastq \
    unpaired${i%?.fastq}2.fastq \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:25 MINLEN:50
done

mkdir ./../results/trimmed

# Then is needed to run Fastqc again in order to check the quality of the files after the cleanning.

for i in ./../results/trimmed*;
do fastqc $i;
done

mv *.zip *.html ./../results/FastQC_Reports/
