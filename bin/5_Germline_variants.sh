# Fot the germline variants calling we are going to continue from the 4_Germline ancestry.sh (line 32).

mkdir ./../results/Germline
Ref="/labs/lgfc/hg19/"
germ="./../results/Germline/" # Directory for germline results 

time gatk GenotypeGVCFs \
    -R $Ref'Homo_sapiens_assembly19.fasta' \
    -V all_gvcf_hg19.g.vcf.gz \
    -O $germ'all_ggvcf_raw_variants.vcf.gz' \
    --allow-old-rms-mapping-quality-annotation-data 2>$germ'GGVCF.log'

time gatk SelectVariants \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V $germ'all_ggvcf_raw_variants.vcf.gz' \
	-select-type SNP \
	-O $germ'all_raw_snps.vcf'

time gatk SelectVariants \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V $germ'all_ggvcf_raw_variants.vcf.gz' \
	-select-type INDEL \
	-O $germ'all_raw_indel.vcf'

time gatk --java-options "-Xmx32g" VariantRecalibrator \
    -R $Ref'Homo_sapiens_assembly19.fasta' \
    -V $germ'all_raw_snps.vcf' \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $Ref'hapmap_3.3.b37.vcf.gz' \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 $Ref'1000G_omni2.5.b37.vcf.gz' \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $Ref'1000G_phase1.snps.high_confidence.b37.vcf.gz' \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Ref'Homo_sapiens_assembly19.dbsnp138.vcf' \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O $germ'all_snp_hg19.recal' \
    --tranches-file $germ'all_snp_hg19.tranches' \
    --rscript-file $germ'all_snp_hg19.plots.R' 2>$germ'Vrecal_snp.log'

time gatk --java-options "-Xmx32g" ApplyVQSR \
    -R $Ref'Homo_sapiens_assembly19.fasta' \
    -V $germ'all_raw_snps.vcf' \
    --truth-sensitivity-filter-level 99.5 \
    --recal-file $germ'all_snp_hg19.recal' \
    --tranches-file $germ'all_snp_hg19.tranches'  \
    -mode SNP \
    -O $germ'all_filt_snps.vcf' 2>$anc'VQSR_indel.log'

time gatk --java-options "-Xmx32g" VariantRecalibrator \
    -R $Ref'Homo_sapiens_assembly19.fasta' \
    -V $germ'all_raw_indel.vcf' \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 $Ref'Mills_and_1000G_gold_standard.indels.b37.sites.vcf' \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Ref'Homo_sapiens_assembly19.dbsnp138.vcf' \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff\
    -mode INDEL \
    -O $germ'all_indel_hg19.recal' \
    --tranches-file $germ'all_indel_hg19.tranches' \
    --rscript-file $germ'all_indel_hg19.plots.R' 2>$germ'Vrecal_indel.log'

time gatk --java-options "-Xmx32g" ApplyVQSR \
    -R $Ref'Homo_sapiens_assembly19.fasta' \
    -V $germ'all_raw_indel.vcf' \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file $germ'all_indel_hg19.recal' \
    --tranches-file $germ'all_indel_hg19.tranches'  \
    -mode INDEL \
    -O $germ'all_filt_indel.vcf' 2>$germ'VQSR_indel.log'

time gatk --java-options '-Xmx64g' VariantFiltration \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V $germ'all_filt_snps.vcf' \
	--filter-expression "QD < 2.0" --filter-name "QD2" \
	--filter-expression "FS > 60.0" --filter-name "FS60" \
	--filter-expression "MQ < 30.0" --filter-name "MQ30" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQRS12.5" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS8" \
	-O $germ'all_filt2_snps.vcf' 2>$germ'VF_snps.log'

time gatk --java-options '-Xmx64g' VariantFiltration \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V $germ'all_filt_indel.vcf' \
	--filter-expression "QD < 2.0" --filter-name "QD2" \
	--filter-expression "FS > 200.0" --filter-name "FS200" \
	--filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS20" \
	-O $germ'all_filt2_indel.vcf' 2>$germ'VF_indels.log'

time gatk --java-options '-Xmx64g' MergeVcfs \
	-I $germ'all_filt2_indel.vcf' \
	-I $germ'all_filt2_snps.vcf'\
    -O $germ'all_filt2_S_I.vcf' 2>$germ'Merge.log'

time gatk --java-options '-Xmx64g' SelectVariants \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V $germ'all_filt2_S_I.vcf' \
	-O $germ'all_germinal_filtered.PASS.vcf' \
    --exclude-filtered 2>$germ'SV_PASS.log'

time gatk Funcotator \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V $germ'all_germinal_filtered.PASS.vcf' \
    -O $germ'all_germinal_annotated_filtered.PASS.vcf' \
	--output-file-format MAF \
	--data-sources-path /labs/lgfc/funcotator_dataSources.v1.7.20200521s \
	--ref-version hg19 2>$germ'germ_funcotator.log'

# VEP can also be used in command line or web tool: https://grch37.ensembl.org/Homo_sapiens/Tools/VEP 