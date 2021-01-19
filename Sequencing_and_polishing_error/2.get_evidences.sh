cat list|grep -v '^#' |while read a
do 
	cd $a
	echo $a
	echo "#python3 ../../bin/select_pileup_from_comparedTab.py $a.pileup $a.compared.tab >$a.10x.compared.pileup
#/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Assess_data_availablity/filter_maternal_reads/Alignmentout_Sum.pl -in $a.10x.compared.pileup -type pileup -ref mat_$a.fa -iTools /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/iTools -outpre $a.10x.compared
#python3 ../../bin/select_pileup_from_comparedTab.py mat2mat.pb.cor.$a.pileup $a.compared.tab >$a.mat2mat.pb.cor.compared.pileup
#python3 ../../bin/select_pileup_from_comparedTab.py mat2mat.pb.raw.$a.pileup $a.compared.tab >$a.mat2mat.pb.raw.compared.pileup
# extract pileup file by coords
#python3 ../../bin/select_Patpileup_from_coords.py pat2pat.pb.cor.$a.pileup $a.rqpos.coords >$a.pat2pat.pb.cor.compared.pileup
#python3 ../../bin/select_Patpileup_from_coords.py pat2pat.pb.raw.$a.pileup $a.rqpos.coords >$a.pat2pat.pb.raw.compared.pileup
#/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Assess_data_availablity/filter_maternal_reads/Alignmentout_Sum.pl -in $a.mat2mat.pb.cor.compared.pileup -type pileup -ref mat_$a.fa -iTools /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/iTools -outpre $a.mat2mat.pb.cor.compared
#/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Assess_data_availablity/filter_maternal_reads/Alignmentout_Sum.pl -in $a.mat2mat.pb.raw.compared.pileup -type pileup -ref mat_$a.fa -iTools /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/iTools -outpre $a.mat2mat.pb.raw.compared
#/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Assess_data_availablity/filter_maternal_reads/Alignmentout_Sum.pl -in $a.pat2pat.pb.cor.compared.pileup -type pileup -ref pat_$a.fa -iTools /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/iTools -outpre $a.pat2pat.pb.cor.compared
#/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Assess_data_availablity/filter_maternal_reads/Alignmentout_Sum.pl -in $a.pat2pat.pb.raw.compared.pileup -type pileup -ref pat_$a.fa -iTools /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/iTools -outpre $a.pat2pat.pb.raw.compared 
#python3 ../../bin/join_evidence.py $a.compared.tab $a.rqpos.coords $a.10x.compared.base.xls $a.mat2mat.pb.raw.compared.base.xls $a.mat2mat.pb.cor.compared.base.xls $a.pat2pat.pb.raw.compared.base.xls $a.pat2pat.pb.cor.compared.base.xls > $a.allinfo.xls
python3 ../../bin/join_evidence1.py $a.compared.tab $a.rqpos.coords $a.10x.compared.base.xls $a.mat2mat.pb.cor.compared.base.xls $a.pat2pat.pb.cor.compared.base.xls > $a.allinfo.xls" >$a.work2.sh
	#sel-qsub evo 0.1g 1 $a.work2.sh
	sh $a.work2.sh
	cd ..
done
