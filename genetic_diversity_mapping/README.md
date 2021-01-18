### Maternal assembly as the reference

#### bwa mem alignment

sample.list contains the sequencing lib IDs:

```text
S1
S2
S3
S4
```


run baw mem alignment

```
ref="data/mat.genome.fa"
gatk="software/gatk-4.1.4.1/gatk"
samtools="software/samtools-1.8/bin/samtools"
bwa="software/bwa"
cat ../sample.list|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	outdir=`pwd`
	echo "$bwa mem -t 10  $ref $data/${a}_1.fq.gz $data/${a}_2.fq.gz | $samtools view -Sb - > $outdir/${a}.bam && echo \"** bwa mapping done **\"
$samtools sort -@ 10 -m 4G -O bam -o $outdir/${a}.sorted.bam $outdir/${a}.bam && echo \"** BAM sort done\"
$gatk MarkDuplicates -I $outdir/${a}.sorted.bam -O $outdir/${a}.sorted.markdup.bam -M $outdir/${a}.sorted.markdup_metrics.txt && echo \"** markdup done **\"
$samtools index $outdir/${a}.sorted.markdup.bam && echo \"** index done **\" " >$a.bwa.sh
	sel-qsub evo 4g 10 $a.bwa.sh
	cd ..
done
```

#### combine sub bam files together, and prepare for variation calling

```
$samtools addreplacerg -@ 4 -r '@RG\tID:foo\tPL:illumina\tSM:OFFSP' -o out.sorted.addRG.bam maternal.merge.sorted.bam
$samtools index out.sorted.addRG.bam
```

#### prepare index for reference
```
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
```

#### run GATK to call variations

##### generating gvcf file

```
time $gatk HaplotypeCaller \
	-R $ref \
	--emit-ref-confidence GVCF \
	-I out.sorted.addRG.bam \
	--native-pair-hmm-threads 12 \
	-O out.g.vcf && echo "** gvcf done **"
```

##### 2 variation calling

```
time $gatk GenotypeGVCFs \
	-R $ref \
	-V out.g.vcf \
	-O out.vcf && echo "** vcf done **"
```
