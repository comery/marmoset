data="/hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/Chrs"
cat list|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	echo "#/share/app/MUMmer-3.23/nucmer --maxmatch -c 500 -b 500 -l 100 $data/$a/mat_$a.fa $data/$a/pat_$a.fa -p $a
#/share/app/MUMmer-3.23/delta-filter -m -i 90 -l 100 $a.delta > $a.delta.filt
#/share/app/MUMmer-3.23/show-coords -THrd $a.delta.filt > $a.delta.filt.coords
source /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Alignment_nomask/syri/conda_zhouyang_init
conda activate syri_env
export PATH=/share/app/MUMmer-3.23/:\$PATH
/hwfssz1/ST_DIVERSITY/PUB/USER/zhouyang/software/syri-1.0/syri/bin/syri -c $a.delta.filt.coords -r $data/$a/mat_$a.fa -q $data/$a/pat_$a.fa -d $a.delta.filt" > $a.syri.sh
	sel-qsub evo 2g 1 $a.syri.sh
	cd ..
done

cat list|while read a
do 
	cd $a
	python3 /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/SVs/Check_syri_nomask/subdivide_SV.py syri.out $a
	cd ..
done
