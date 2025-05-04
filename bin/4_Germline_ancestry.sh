# Ancestry determination.

Res="/labs/lgfc/Results/"
Ref="/labs/lgfc/hg19/"

for i in $Res*'_ready.bam'
do
	time gatk --java-options "-Xmx4g" HaplotypeCaller \
	    -R $Ref'Homo_sapiens_assembly19.fasta'\
	    -I $i \
	    -L /labs/lgfc/truseq_nochr.bed \
	    -O ${i%_ready.bam}'.g.vcf.gz' \
	    -ERC GVCF 2>${i%_ready.bam}_HC.log
done

# We need to include all the samples (P01, P02, P03,..., Pn)
time gatk4 CombineGVCFs \
    -R $Ref'Homo_sapiens_assembly19.fasta' \
	--variant $Res'P01.g.vcf.gz' \
	--variant $Res'P02.g.vcf.gz' \
	--variant $Res'P03.g.vcf.gz' \
	--variant $Res'P04.g.vcf.gz' \
	-O all_gvcf_hg19.g.vcf.gz \
	--dbsnp $Ref'Homo_sapiens_assembly19.dbsnp138.vcf' 2>combineGVCFs.log

time gatk4 --java-options "-Xmx32g" GenotypeGVCFs \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V all_gvcf_hg19.g.vcf.gz \
	-O all_ggvcf_hg19.vcf.gz \
	--allow-old-rms-mapping-quality-annotation-data True \
	--dbsnp $Ref'Homo_sapiens_assembly19.dbsnp138.vcf' 2>GGVCF.log

## FROM THIS POINT IS POSSIBLE TO CONTINUE TO GERMLINE VARIANTS DETERMINATION (5_Germline_variants.sh)

time gatk4 --java-options "-Xmx32g" VariantRecalibrator \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V all_ggvcf_hg19.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $Ref'hapmap_3.3.b37.vcf.gz' \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 $Ref'1000G_omni2.5.b37.vcf.gz' \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 $Ref'1000G_phase1.snps.high_confidence.b37.vcf.gz' \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $Ref'Homo_sapiens_assembly19.dbsnp138.vcf' \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-O all_hg19.recal \
	--tranches-file all_hg19.tranches \
	--rscript-file all_hg19.plots.R 2>Vrecal.log


time gatk4 --java-options "-Xmx32g" ApplyVQSR \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-V all_ggvcf_hg19.vcf.gz \
	--truth-sensitivity-filter-level 99.0 \
	--recal-file all_hg19.recal \
	--tranches-file all_hg19.tranches \
	-mode SNP \
	-O all_VQSR_hg19.vcf.gz 2>VQSR.log


vcftools \
	--gzvcf all_VQSR_hg19.vcf.gz \
	--remove-filtered-all \
	--max-missing 0.98 \
	--recode \
	--out all_PASS_hg19.vcf.gz 2>vcftools_PASS.log


plink \
	--vcf all_PASS_hg19.vcf.gz.recode.vcf \
	--double-id \
	--keep-allele-order \
	--make-bed \
	--out all_PASS_hg19


plink \
	--bfile AFR-EUR-AME_autosomes_seq.geno98 \
	--bmerge  all_PASS_hg19.bed all_PASS_hg19.bim all_PASS_hg19.fam \
	--geno 0.05 \
	--make-bed \
	--out Pops_LUAD_snps-only

plink \
    --bfile all_PASS_hg19 \
    --exclude Pops_LUAD_snps-only-merge.missnp\
    --make-bed \
    --out LUAD_exclude

plink \
    --allow-no-sex \
	--bfile $Ref'AFR-EUR-AME_autosomes_seq.geno98' \
    --bmerge  LUAD_exclude.bed LUAD_exclude.bim LUAD_exclude.fam \
	--geno 0.05 \
    --make-bed \
    --out REF_POPS_LUAD

admixture REF_POPS_LUAD.bed 3 > REF_POPS_LUAD_admixture.log
mv REF_POPS_LUAD.3.* $Res

paste REF_POPS_LUAD.fam REF_POPS_LUAD.3.Q > LUAD_admixture.3.Q.ids.txt