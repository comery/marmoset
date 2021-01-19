data='data/SRP041711'
cat  SraRunInfo.txt|grep -v '^#' |while read a b
do
	[ -d $b ] || mkdir $b
	cd $b
	ln -s $data/$a/*.sra 
	for i in `ls *.sra`
	do
		echo "/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/NGS_reads_tools/sratoolkit.2.9.2/bin/fastq-dump --split-3 $i" 
	done >sra2fq.sh
	# run sra2fq.sh
	cd ..
done


# fq2clean
cat list|while read a
do
	cd $a
	ls *.sra|sed 's/.sra//g'|while read b
	do
		echo "/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/NGS_reads_tools/fastp -i ${b}_1.fastq -I ${b}_2.fastq -o ${b}_1.clean.fq.gz -O ${b}_2.clean.fq.gz && rm ${b}_1.fastq ${b}_2.fastq "
	done > fq2clean.sh
	# run fq2clean.sh
	cd ..
done

# bwa mem mapping
ref="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/0.assembly/curated/mat.genome.curated.fa"
cat list|while read a
do
	cd $a
	ls *.sra|sed 's/.sra//g'|while read b
	do
		echo "bwa mem -t 6  $ref ${b}_1.clean.fq.gz ${b}_2.clean.fq.gz | samtools view -q 10 -Sb - > ${a}.bam && echo \"** $b mapping done *"
	done > bwa_mem.sh
	echo "bam=\`find ./ -name \"*.bam\"\` " >> bwa_mem.sh
	echo "samtools merge -@ 6 $a.merge.bam \$bam && echo \"** merge done **\" " >>bwa_mem.sh
	echo "samtools sort -@ 6 -m 2G -O bam -o ${a}.sorted.bam \$bam && echo \"** BAM sort done **\"" >> bwa_mem.sh
	echo "samtools index ${a}.sorted.bam && echo \"** index done **\"" >> bwa_mem.sh
	# run bwa_mem.sh
	cd ..
done


# call snp
ref="data/mat.genome.fa"
gatk="software/gatk-4.1.4.1/gatk"
samtools="software/samtools-1.8/bin/samtools"
bcftools="software/bcftools-1.8/bin/bcftools"
cat list|while read a
do
	cd $a
	echo "#$samtools view  -h  -F 12 -F 256 $a.sorted.bam |grep -v 'XA:Z' | grep -v 'SA:Z' |$samtools view -b - >$a.sorted.uniq.bam
#$samtools addreplacerg  -r '@RG\tID:$a\tPL:illumina\tSM:$a' -o $a.sorted.uniq.addRG.bam $a.sorted.uniq.bam
$samtools index $a.sorted.uniq.addRG.bam
$gatk HaplotypeCaller \\
	-R $ref \\
	--emit-ref-confidence GVCF \\
	-I $a.sorted.uniq.addRG.bam \\
	--native-pair-hmm-threads 10 \\
	-O $a.g.vcf && echo \"** gvcf done **\"
$gatk GenotypeGVCFs \\
	-R $ref \\
	-V $a.g.vcf \\
	-O $a.vcf && echo \"** vcf done **\" " >$a.gatk.sh
	#run $a.gatk.sh
	cd ..
done

# filter snp
gatk="software/gatk-4.1.4.1/gatk"
cat list|while read a
do
	cd $a
	$gatk SelectVariants \
		-select-type SNP \
		-V $a.vcf \
		-O $a.snp.vcf.gz

	$gatk VariantFiltration \
		-V $a.snp.vcf.gz \
		--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
		--filter-name "Filter" \
		-O $a.snp.filter.vcf.gz
	cd ..
done


# keep chr1-chr22 and claspe filtered SNP
cat list|while read a
do
	cd $a
	gzip -dc $a.snp.filter.vcf.gz|grep -v -e 'scaffold' -e 'Filter' |sed 's/mat_//g' > $a.snp.good.vcf
	bgzip $a.snp.good.vcf
	tabix $a.snp.good.vcf.gz
	/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Evolution/plink-v1.90b5.1/plink --vcf $a.snp.good.vcf --maf 0.05 --recode --out $a.plink
	cd ..
done


##draw plot for each individual after obtain the depth
# the input.txt contains 5 columns "Sample\tChromosome\tLocation\tDepth\tLogNum"
library(ggplot2)
a<-read.table("input.txt",sep="\t",header=F)
b<-log(a$Depth)
c<-cbind(a,b)
colnames(c)<-c("Sample","Chromosome","Location","Depth","LogNum")
pal<-c('#E9B7ED','#f38181','#fce38a','#00b8a9')
c$Chromosome<- factor(c$Chromosome,levels=c("Chr2","Chr3","Chr4","Chr13"))
ggplot(data=c,aes(x=Location,y=LogNum,fill=Sample,color=Sample,group=Sample))+geom_line()+facet_grid(Chromosome~.)+xlab("Location")+ylab("Log(Variants_Number)")+scale_size_manual(values=c(0.1, 0.1,0.1,0.1))+ scale_colour_manual(values = pal)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_set(theme_bw())
dev.off()


# Roh analysis using all g.vcf
#combine gvcfs for 11 individuals
list_of_gvcfs="gvcf_11samples.lis"
##"gvcf_11samples.lis": the list file with 11 IDs.out.g.vcf
gvcf_argument=""
while read -r line; do
gvcf_argument=$gvcf_argument" -v $line"
done < "$list_of_gvcfs"

$sentieon driver -r $ref   --algo GVCFtyper  \
 --call_conf  30 --emit_conf 30  \
 --emit_mode variant  \
 --max_alt_alleles 2 \
 --annotation 'AC,AF,AN,BaseQRankSum,QD,FS,SOR,MQ,MQRankSum,ReadPosRankSum' \
 $gvcf_argument joint.vcf.gz

#SNP calling for 11 samples
time gatk-4.0.5.2/gatk SelectVariants \
    -select-type SNP \
    -V joint.vcf.gz \
    -O snp.vcf.gz
time gatk-4.0.5.2/gatk SelectVariants \
    -select-type INDEL \
    -V joint.vcf.gz \
    -O indel.vcf.gz
time gatk-4.0.5.2/gatk VariantFiltration \
    -V indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "HARD_fail" \
    -O indel.filter.vcf.gz

time gatk-4.0.5.2/gatk VariantFiltration \
    -V snp.vcf.gz \
    --filter-expression "AN < 22" --filter-name "AN_FAIL" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL_FAIL" \
    --filter-expression "DP < 50.0 || DP > 1500.0" --filter-name "DP_FAIL" \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "HARD_FAIL" \
    --mask indel.filter.vcf.gz --mask-extension 3 --mask-name aroundIndel \
    --cluster-size 3 --cluster-window-size 10 \
    -O snp.filter.vcf.gz
time gatk-4.0.5.2/gatk SelectVariants \
    -V snp.filter.vcf.gz \
    --exclude-filtered \
    -O snp.pass.vcf.gz

vcftools --max-missing 1  --recode-INFO-all  --recode --out final --gzvcf snp.pass.vcf.gz   --min-alleles 2 --max-alleles 2  --remove-filtered-all   --stdout | bgzip -c > final.snp.vcf.gz
vcftools --FILTER-summary  --out check --gzvcf final.snp.vcf.gz

# offspring's snp was generated by mummer alignment, so here you need change its format into vcf.

#combine all 12 samples
bcftools-1.8/bin/bcftools merge --output combine.vcf.gz --output-type z --missing-to-ref  --threads 2 final.snp.vcf.gz  offspring.vcf.gz
#keep only normal chromosomes
gzip -dc combine.vcf.gz  |grep -v ChrX | grep -v scaffold |bgzip -c >combine_noX.vcf.gz
#maf & hwe
vcftools --max-missing 1  --maf 0.2  --hwe 0.001  --recode-INFO-all  --recode  --gzvcf  combine_noX.vcf.gz  --stdout |bgzip -c > use.vcf.gz
#ROH
software/bcftools-1.8/bcftools roh  -O r --AF-dflt 0.4 -o roh  use.vcf.gz
awk '$7>=10' roh>roh.fil.results

