# Fot the somatic variants calling of one tumor sample.
# We need to have our reference genome.

mkdir ./../results/Somatic
Res=" ./../results/Somatic/"
Ref="/labs/lgfc/hg19/"

for i in $Res*'_T_ready.bam'
do
name="${i##*/}"
nombre="${name%%'_'*}"
nombreN="${nombre%%'T'}""N"
time gatk  --java-options "-Xmx16G" Mutect2 \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	-I $i \
	-I ${i%_T_ready.bam}'_N.bam' \
	-O ${i%_T_ready.bam}.vcf.gz \
	-tumor $nombre \
	-normal $nombreN \
	--germline-resource $Ref'af-only-gnomad.raw.sites.grch37.vcf.gz' \
	--panel-of-normals $Ref'PON_hg19_x3.vcf.gz' \
	--f1r2-tar-gz ${i%_T_ready.bam}'_f1r2.tar.gz' \
	-mbq 20 2>${i%_T_ready.bam}'_m2.log' 	

time gatk LearnReadOrientationModel \
	-I ${i%_T_ready.bam}'_f1r2.tar.gz' \
	-O ${i%_T_ready.bam}'_rom.tar.gz'  2>${i%_T_ready.bam}_lom.log
	
time gatk GetPileupSummaries -I ${i%_T_ready.bam}'N.bam' \
	-V $Ref'ExAC.r1.sites.vep.vcf.gz' \
	-L $Ref'ExAC.r1.sites.vep.vcf.gz' \
	-O ${i%_T_ready.bam}'N_ps.table' 2>${i%_T_ready.bam}_gpsN.log
	
time gatk GetPileupSummaries -I $i \
	-V $Ref'ExAC.r1.sites.vep.vcf.gz' \
	-L $Ref'ExAC.r1.sites.vep.vcf.gz' \
	-O ${i%_T_ready.bam}'_ps.table' 2>${i%_T_ready.bam}_gpsT.log	

time gatk CalculateContamination \
	-I ${i%_T_ready.bam}'_ps.table' \
	--matched ${i%_T_ready.bam}'N_ps.table' \
   	--tumor-segmentation ${i%_T_ready.bam}'_segments.table' \
	-O ${i%_T_ready.bam}'_cc.table'  2>${i%_T_ready.bam}_cc.log
	
time gatk FilterMutectCalls \
	-V ${i%_T_ready.bam}.vcf.gz \
	-R $Ref'Homo_sapiens_assembly19.fasta' \
	--contamination-table ${i%_T_ready.bam}'_cc.table' \
	--ob-priors ${i%_T_ready.bam}'_rom.tar.gz' \
	--tumor-segmentation ${i%_T_ready.bam}'_segments.table'  \
	-O ${i%_T_ready.bam}'_filt.vcf.gz'  2>${i%_T_ready.bam}_f2.log

bcftools view \
    -f PASS ${i%_T_ready.bam}'_filt.vcf.gz'  > ${i%_T_ready.bam}'_PASS.vcf'
done

for i in $Res*'PASS.vcf'
do
      time gatk4 Funcotator -R $Ref'Homo_sapiens_assembly19.fasta' \
	       -V $i \
	       -O ${i%.vcf}'_annotated.maf' \
	       --output-file-format MAF \
	       --data-sources-path funcotator_dataSources.v1.7.20200521s \
	       --ref-version hg19 2>${i%_PASS.vcf}'_funcotator.log'
done

# VEP option for annotation

for i in $Res*'_PASS.vcf'
do
    vep --cache --cache_version 110 --dir /labs/lgfc/hg19 \
	    --port 3337 \
	    -i \
	    -o ${i%_PASS.vcf}'_PASS_vep.vcf' \
	    --everything 2>${i%_filt.vcf.gz}_vep.log
done