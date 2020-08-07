data='/hwfssz1/pub/database/ftp.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP041/SRP041711'
cat  SraRunInfo.txt|grep -v '^#' |while read a b
do
	[ -d $b ] || mkdir $b
	cd $b
	#ln -s $data/$a/*.sra 
	#for i in `ls *.sra`
	#do
	#	echo "/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/NGS_reads_tools/sratoolkit.2.9.2/bin/fastq-dump --split-3 $i" 
	#done >sra2fq.sh
	cd ..
done
#sra2fq
#cat list|while read a;do cd $a; sel-qsub evo 0.1g 1 sra2fq.sh;cd ..;done

# fq2clean
cat list|while read a
do
	cd $a
	ls *.sra|sed 's/.sra//g'|while read b
	do
		echo "/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/NGS_reads_tools/fastp -i ${b}_1.fastq -I ${b}_2.fastq -o ${b}_1.clean.fq.gz -O ${b}_2.clean.fq.gz && rm ${b}_1.fastq ${b}_2.fastq "
	done > fq2clean.sh
	sel-qsub evo 2g 2 fq2clean.sh
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
	# sel-qsub evo 2g 6 bwa_mem.sh
	cd ..
done


# call snp
ref="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/0.assembly/curated/mat.genome.curated.fa"
gatk="/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/gatk-4.1.4.1/gatk"
samtools="/hwfssz5/ST_DIVERSITY/B10K/PUB/local/basic_bio-tool/samtools-1.8/bin/samtools"
bcftools="/hwfssz5/ST_DIVERSITY/B10K/PUB/local/basic_bio-tool/bcftools-1.8/bin/bcftools"
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
	#sel-qsub evo 6g 10 $a.gatk.sh
	cd ..
done

# filter snp
gatk="/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/gatk-4.1.4.1/gatk"
cat list|while read a
do
	cd $a
	# 使用SelectVariants，选出SNP
	$gatk SelectVariants \
		-select-type SNP \
		-V $a.vcf \
		-O $a.snp.vcf.gz

	# 为SNP作硬过滤
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


# merge all vcf
bcftools merge -f PASS -O z SRS602594/SRS602594.snp.filter.vcf.gz offspring/all.diff.snp.artifical.vcf.gz
