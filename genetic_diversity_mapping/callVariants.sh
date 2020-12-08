# maternal assembly will be the reference
ref="data/mat.genome.fa"
gatk="software/gatk-4.1.4.1/gatk"
samtools="software/samtools-1.8/bin/samtools"
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
#1 generating gvcf file
time $gatk HaplotypeCaller \
	-R $ref \
	--emit-ref-confidence GVCF \
	-I out.sorted.addRG.bam \
	--native-pair-hmm-threads 12 \
	-O out.g.vcf && echo "** gvcf done **"

#2 variation calling
time $gatk GenotypeGVCFs \
	-R $ref \
	-V out.g.vcf \
	-O out.vcf && echo "** vcf done **"
