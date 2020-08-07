cat list |while read a
do 
	[ -d $a  ] || mkdir $a
	cd $a
	#ln -s ../../Mummer/$a/$a.diff.snp .
	#ln -s ../../Mummer/$a/$a.diff.1coords .
	#ln -s ../../GATK_10X/$a/$a.pileup .
	#ln -s ../../GATK_10X/$a/$a.sambcf.vcf.gz
	#ln -s ../../GATK_10X/$a/$a.snp.vcf.gz
	#ln -s ../../../Alignment_nomask/Chrs/$a/mat_$a.fa
	#ln -s ../../../Alignment_nomask/Chrs/$a/pat_$a.fa
	#ln -s /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/tmp/marmoset/Mapping_Pacbio/pat2pat/pat_Chrs/pat2pat.pb.cor.$a.pileup .
	#ln -s /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/tmp/marmoset/Mapping_Pacbio/mat2mat/mat_Chrs/mat2mat.pb.cor.$a.pileup .
	#rm mat2mat.pb.raw.$a.pileup pat2pat.pb.raw.$a.pileup
	#ln -s /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/tmp/marmoset/raw_pacbio/mat/mat_Chrs/mat2mat.pb.raw.$a.pileup .
	#ln -s /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/tmp/marmoset/raw_pacbio/pat/pat_Chrs/pat2pat.pb.raw.$a.pileup .
	ln -s /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Mummer/$a/$a.rqpos.coords .
	cd ..
done
