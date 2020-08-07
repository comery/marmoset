# split genome by chrs
#python3 split_chrs.py /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/0.assembly/Sanger_20191216/mCalJac1.maternal.chrs.fa /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/0.assembly/Sanger_20191216/mCalJac1.paternal.chrs.fa

# make up vcf file and statistic snp density
cat list|while read a
do
	cd $a
	#echo "#/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/All_kinds_align/MUMmer-3.23/nucmer -maxmatch -l 100 -c 500 -p $a mat_$a.fa pat_$a.fa
#/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/All_kinds_align/MUMmer-3.23/dnadiff -d $a.delta -p $a.diff
#mv $a.diff.snps $a.diff.var
#python3 ../split_snp_indel.py $a.diff.var $a.diff" >$a.work.sh
	#sel-qsub evo 0.5g 1 $a.work.sh
	
	# change indel format for easy reading
	#python3 indel_statistic.py $a.diff.indel $a
	
	# split small and longer indels
	#awk '$6<=50' $a.diff.indel.report.txt >$a.mummer.sm.indel
	#awk '$6>50' $a.diff.indel.report.txt >$a.mummer.lg.indel
	
	# small indel supported by 10x at least 3 reads
	#python3 ../filter_indel_by_depth.py $a.mummer.sm.indel /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/SNP_INDEL/GATK_10X/$a/$a.10x.base.xls  $a.mummer.sm.10xsupport.indel 1>sm.indel.check.log
	
	# large indel supported by Pacbio long reads at least 3 reads
	#echo "python3 ../filter_indel_by_depth.py $a.mummer.lg.indel /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/tmp/marmoset/corrected_pacbio/pat2mat/$a/$a.pat2mat.cor.base.xls  $a.mummer.lg.Pbsupport.indel 1>lg.indel.check.log " > $a.check.lg.sh
	#sel-qsub evo 200m 1 $a.check.lg.sh
	
	## calculate snp density
	### sliding window = 100k, step = 20k
	python3 ../calSNPdensity.py $a.diff.snp ../chr.curated.len 100000 20000 $a >$a.100k_20k.snp.density

	### sliding window = 100k, step = 0
	python3 ../calSNPdensity.py $a.diff.snp ../chr.curated.len 100000 100000 $a >$a.100k_100k.snp.density
	cd ..
done

