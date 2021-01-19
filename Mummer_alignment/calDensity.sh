cat list|while read a
do 
	cd $a
	# window length = 100k, step length = 20k (overlap = 80k)
	python3 ../calSNPdensity.py $a.diff.snp ../chr.curated.len 100000 20000 $a > $a.100k_20k.snp.density
	# window length = 100k, step length = 100k (no overlap)
	python3 ../calSNPdensity.py $a.diff.snp ../chr.curated.len 100000 100000 $a > $a.100k_100k.snp.density
	python3 ../calINDELdensity.py $a.mummer.sm.indel ../chr.curated.len 100000 20000 $a > $a.100k_20k.indel.density
	python3 ../calINDELdensity.py $a.mummer.sm.10xsupport.indel ../chr.curated.len 100000 20000 $a > $a.100k_20k.10xsupport.indel.density
	python3 ../calINDELdensity.py $a.mummer.sm.indel ../chr.curated.len 100000 100000 $a > $a.100k_100k.indel.density
	python3 ../calINDELdensity.py $a.mummer.sm.10xsupport.indel ../chr.curated.len 100000 100000 $a > $a.100k_100k.10xsupport.indel.density

	cd ..
done
