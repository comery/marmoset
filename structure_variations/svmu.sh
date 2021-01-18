dir="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask"
snp_mode="l"
cat ../list|grep -v '^#'|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	#echo "/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/svmu-0.4-alpha/svmu $dir/Mummer_syri/$a/$a.delta.filt $dir/Chrs/$a/mat_$a.fa $dir/Chrs/$a/pat_$a.fa $snp_mode $dir/Lastz/$a/$a.lastz.txt $a && echo done" > $a.svmu.sh
	#sel-qsub evo 2g 1 $a.svmu.sh
	# subtype indels larger than 200 bp defined as  large indels; (< 200 bp indels are more likely to alignment boundary issues)
	grep -v 'REF_CHROM' sv.$a.txt |awk '$4=="INS"' >$a.svmu.lg.ins.txt
	grep -v 'REF_CHROM' sv.$a.txt |awk '$4=="DEL"' >$a.svmu.lg.del.txt
	python3 ../filter_indel_len.py $a.svmu.lg.ins.txt > $a.svmu.lg.ins.filter.txt
	python3 ../filter_indel_len.py $a.svmu.lg.del.txt > $a.svmu.lg.del.filter.txt
	grep -v -e 'REF_CHROM' -e 'INS' -e 'DEL' -e 'CNV' sv.$a.txt >$a.svmu.nonINDEL.txt
	cd ..
done


cat ../list |while read a;do cd $a; cat $a.svmu.lg.ins.filter.txt $a.svmu.lg.del.filter.txt > $a.svmu.lg.indel.filter.txt; python3 ../calINDELdensity.py $a.svmu.lg.indel.filter.txt /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Mummer/chr.curated.len 1000000 1000000 $a >$a.indel.1m.density;cd ..;done
