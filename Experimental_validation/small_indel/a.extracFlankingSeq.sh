indel="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Mummer"
cat ../list|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	[ -e $a.mummer.sm.10xsupport.indel  ] || ln -s $indel/$a/$a.mummer.sm.10xsupport.indel
	[ -e $a.alns ] || ln -s /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Mummer/$a/$a.alns
	[ -e $a.repeat.bed ] || ln -s ../../1.select_snp_for_validation/$a/$a.repeat.bed
	echo "#python3 ../indelPos2bed.py $a.mummer.sm.10xsupport.indel >indel.target.bed
# extend target position upto template length, template_length = 500
#awk '{print \"mat_\"\$1\"\\t\"\$2-250\"\\t\"\$3+250}' indel.target.bed >indel.temp.bed
#bedtools intersect -a indel.temp.bed -b $a.repeat.bed -loj >indel.temp.repeat.overlap.bed
#grep '\\.' indel.temp.repeat.overlap.bed|awk '{print \$1\"\\t\"\$2+250}' >indel.temp.outrepeat.pos
#python3 ../extractPos_from_bed.py indel.temp.outrepeat.pos $a.mummer.sm.10xsupport.indel >all.potenial.pos
python3 ../extractTemplate_from_MummerAlns.py -in $a.alns -pos all.potenial.pos -out all.potenial.template.fa -len 250 -extra_snp_allow 1" >$a.extractTemp.sh
	sel-qsub evo 0.1g 1 $a.extractTemp.sh
	#sh $a.extractTemp.sh
	cd ..
done
