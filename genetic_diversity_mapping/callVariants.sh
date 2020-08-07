ref="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/0.assembly/Sanger_20191216/2.chrs/mat.genome.fa"
gatk="/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/gatk-4.1.4.1/gatk"
samtools="/hwfssz5/ST_DIVERSITY/B10K/PUB/local/basic_bio-tool/samtools-1.8/bin/samtools"
#0 check ref's index fai
wkdir=`pwd`
if [ -f $ref.fai ]
then
	echo "fasta index fai of reference is fine!"
else
	cd ${ref%/*}
	$samtools faidx $ref
	cd $wkdir
fi

if [ -f ${ref%%.*}.dict ]
then
	echo "fasta index dict of reference is fine!"
else
	cd ${ref%/*}
	$gatk CreateSequenceDictionary -R $ref
	cd $wkdir
fi
$samtools addreplacerg -@ 4 -r '@RG\tID:foo\tPL:illumina\tSM:OFFSP' -o out.sorted.addRG.bam maternal.merge.sorted.bam
$samtools index out.sorted.addRG.bam
#1 生成中间文件gvcf
time $gatk HaplotypeCaller \
	-R $ref \
	--emit-ref-confidence GVCF \
	-I out.sorted.addRG.bam \
	--native-pair-hmm-threads 12 \
	-O out.g.vcf && echo "** gvcf done **"

#2 通过gvcf检测变异
time $gatk GenotypeGVCFs \
	-R $ref \
	-V out.g.vcf \
	-O out.vcf && echo "** vcf done **"
