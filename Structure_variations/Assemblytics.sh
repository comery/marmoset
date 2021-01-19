data="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Mummer"
assemblytics="/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/Assemblytics-1.2/scripts/Assemblytics"
cat ../list|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	[ -e $a.delta ] || ln -s $data/$a/$a.delta
	# Assemblytics delta output_prefix unique_length_required min_size max_size
	#$assemblytics $a.delta $a 10000 50 100000 && echo "$a done!"
	#/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/Assemblytics-1.2/scripts/Assemblytics_dotplot.R $a
	python3 /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Genome/Assemblytics-1.2/scripts/Assemblytics_within_alignment.py --delta $a.delta --min 50 --output $a.variants_within_alignments.bed
	cat $a.variants_within_alignments.bed $a.variants_between_alignments.bed|sort -k2n > $a.Assemblytics_structural_variants.bed

	cd ..
done

# filter by length and N rate
ref="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Chrs/$a"
cat ../list|while read a
do
	cd $a
	python3 ../check_assemblytics_variants.py $a.Assemblytics_structural_variants.bed $ref/mat_$a.fa  $ref/pat_$a.fa  >$a.Assemblytics_structural_variants.gc.n.bed
	mv $a.Assemblytics_structural_variants $a.cnv.gc.n.bed
	awk '$4>=10 && $10>=10' $a.cnv.gc.n.bed|awk '$6<0.5 && $12<0.5' >$a.cnv.fil.bed
	cd ..
done

# statistic 
cat ../list |while read a;do cat $a/$a.cnv.fil.bed ;done >all.cnv.fil.bed

cat ../list |while read a;do cd $a; python3 ../calCNVdensity.py $a.cnv.fil.bed /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Mummer/chr.curated.len 1000000 1000000 $a >$a.cnv.1m.density;cd ..;done
