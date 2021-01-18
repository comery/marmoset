### SNP, indel, and SV detection using whole haplotype genome alignment

To call heterozygous sites between the two haploid sequences, independent of the GenomeScope calculation, we first performed a Mummer alignment. Because our assemblies span most repetitive sequences, repeat-masking treatment was not necessary before conducting Mummer alignment. 

```
split_genome_bychrs.py
```

### split genome by chrs

two haplotype assemblies will be splited into different folder by chromosome ID.

```shell
python3 split_chrs.py mCalJac1.maternal.chrs.fa mCalJac1.paternal.chrs.fa
```
Note: the fasta sequence header is formatted like 'mat_Chr1' for maternal, and 'patChr1' for maternal.

### Mummer alignment

```shell
cat list|while read a
do
	cd $a
	echo "software/MUMmer-3.23/nucmer -maxmatch -l 100 -c 500 -p $a mat_$a.fa pat_$a.fa
software/MUMmer-3.23/dnadiff -d $a.delta -p $a.diff
mv $a.diff.snps $a.diff.var
python3 ../split_snp_indel.py $a.diff.var $a.diff" >$a.work.sh
	# run $a.work.sh

	# change indel format for easy reading
	python3 indel_statistic.py $a.diff.indel $a

	# split small and longer indels
	awk '$6<=50' $a.diff.indel.report.txt >$a.mummer.sm.indel
	awk '$6>50' $a.diff.indel.report.txt >$a.mummer.lg.indel

	# small indel supported by 10x at least 3 reads
	# here $a.10x.base.xls is calculated from bam file (10x reads mapping to maternal assembly)
	python3 ../filter_indel_by_depth.py $a.mummer.sm.indel xxx/GATK_10X/$a/$a.10x.base.xls  $a.mummer.sm.10xsupport.indel 1>sm.indel.check.log

	# large indel supported by Pacbio long reads at least 3 reads
	# here $a.pat2mat.cor.base.xls is calculated form bam file (paternal pacbio reads mapping to maternal assembly)
	echo "python3 ../filter_indel_by_depth.py $a.mummer.lg.indel xxx/pat2mat/$a/$a.pat2mat.cor.base.xls  $a.mummer.lg.Pbsupport.indel 1>lg.indel.check.log " > $a.check.lg.sh
	# run $a.check.lg.sh

	## calculate snp density
	### sliding window = 100k, step = 20k
	python3 ../calSNPdensity.py $a.diff.snp ../chr.curated.len 100000 20000 $a >$a.100k_20k.snp.density

	### sliding window = 100k, step = 0
	python3 ../calSNPdensity.py $a.diff.snp ../chr.curated.len 100000 100000 $a >$a.100k_100k.snp.density
	cd ..
done
```