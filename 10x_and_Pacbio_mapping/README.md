### Maternal assembly as the reference

#### bwa mem alignment

sample.list contains the sequencing lib IDs:

```text
S1
S2
S3
S4
```


#### run baw mem alignment using 10x illumina reads

```shell
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

#### run ngmlr using corrected Pacbio reads

```shell
ngmlr="software/All_kinds_align/ngmlr-0.2.7/ngmlr"
pat_genome="xxx/pat.genome.curated.fa"
mat_genome="xxx/mat.genome.curated.fa"
$ngmlr -t 8 -r $pat_genome -q xxx/marmoset/corrected_pacbio/pat/pat.pb.corrected.fa.gz |samtools view -bS - >pat.pb.corrected.bam

$ngmlr -t 8 -r $pat_genome -q xxx/marmoset/corrected_pacbio/mat/mat.pb.corrected.fa.gz |samtools view -bS - >mat.pb.corrected.bam

```

#### run ngmlr using raw Pacbio reads

##### 1. assign raw pacbio reads by sequence ids

```
gzip -dc triocanu/mat/asm.correctedReads.fasta.gz |grep '^>' |sed 's/>//g' >mat.asm.correctedReads.id
gzip -dc triocanu/pat/asm-haplotypePaternal.correctedReads.fasta.gz|grep '^>' |sed 's/>//g' >pat.asm.correctedReads.id
python3 assign_raw_Pacbio.py raw.fasta mat.asm.correctedReads.id pat.asm.correctedReads.id raw_pacbio

```

##### 2. run ngmlr

```shell
ngmlr="software/All_kinds_align/ngmlr-0.2.7/ngmlr"
pat_genome="xxx/pat.genome.curated.fa"
mat_genome="xxx/mat.genome.curated.fa"
$ngmlr -t 8 -r $pat_genome -q xxx/marmoset/raw_pacbio/pat/pat.pb.raw.fa.gz |samtools view -bS - >pat.pb.raw.bam

$ngmlr -t 8 -r $pat_genome -q xxx/marmoset/raw_pacbio/mat/mat.pb.raw.fa.gz |samtools view -bS - >mat.pb.raw.bam

```

