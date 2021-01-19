<meta charset="utf-8" emacsmode="-*- markdown -*-"> <link rel="stylesheet" href="https://casual-effects.com/markdeep/latest/journal.css?">



# SVs



Structural variation was originally defined as insertions, deletions and inversions greater than 50 bp in size.


### detection of structure variation

A series of custom scripts (https://github.com/comery/marmoset) identified and sorted our SNPs, indels in the alignments. We employed svmu (v0.4-alpha), Assemblytics (v1.2) and SyRi (v1.0), to detect SVs from Mummer alignment. After several test rounds, we found that svmu reported more accurate large indels, Assemblytics detected CNVs, particularly tandem repeats, whereas SyRi detected other SVs well. Hence, we employed these three methods and combined the results as confident SVs. For Mummer alignment, we used the settings of “nucmer --maxmatch -c 500 -b 500 -l 100” and following filtering of “delta-filter -m -i 90 -l 100”, which was recommended by SyRi (https://schneebergerlab.github.io/syri/), while we used default parameters for svmu and Assemblytics. 


### run SyRi to detect SVs (INV, INVTR, TRANS)

current directory: xxx/marmoset/2.genetic_divergency/Alignment_nomask/syri

The SyRI (Synteny and Rearrangement Identifier), a method to identify structural as well as sequence differences between two whole-genome assemblies. SyRI expects whole-genome alignments (WGA) as input and starts by searching for differences in the structures of the genomes. Afterwards, SyRI identifies local sequence differences within both the rearranged and the non-rearranged (syntenic) regions. SyRI annotates the coordinates of rearranged regions (i.e., breakpoints on both sides of a rearrangement in both genomes) providing a complete regional annotation of rearrangements. 

~~~~~bash
data="xxx/marmoset/2.genetic_divergency/Alignment_nomask/Chrs"
cat ../chr.list|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	echo "$software/MUMmer-3.23/nucmer --maxmatch -c 500 -b 500 -l 100 $data/$a/mat_$a.fa $data/$a/pat_$a.fa -p $a
$software/MUMmer-3.23/delta-filter -m -i 90 -l 100 $a.delta > $a.delta.filt
$software/MUMmer-3.23/show-coords -THrd $a.delta.filt > $a.delta.filt.coords
conda activate syri_env
export PATH=$software/MUMmer-3.23/:\$PATH
$software/syri-1.0/syri/bin/syri -c $a.delta.filt.coords -r $data/$a/mat_$a.fa -q $data/$a/pat_$a.fa -d $a.delta.filt" > $a.syri.sh
	sel-qsub evo 2g 1 $a.syri.sh
	cd ..
done
~~~~~


### run svmu to detect large INDELs

SVMU (Structural Variants from MUmmer) attempts to identify comprehensive sequence variants via alignment of two contiguous genome assemblies. It combines the strengths of two powerful aligners MUMmer and LASTZ to annotate duplicates, large indels, inversions, small indels, SNPs from whole genome alignments.

1. Obtain the mummer and lastz alignments:

Note: Lastz output is not suit for this study, when I add lastz results into svmu input, there is nothing
coming to result, but when I add an empty lastz file, it worked.

mummer alignment delat file can get from above, but here we need do a filtering step.


2. run svmu:

current directory: xxx/marmoset/2.genetic_divergency/Alignment_nomask/svmu

~~~~~bash
dir="xxx/marmoset/2.genetic_divergency/Alignment_nomask"
snp_mode="l"
cat ../chr.list|grep -v '^#'|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	echo "$software/svmu-0.4-alpha/svmu $dir/syri/$a/$a.delta.filt $dir/Chrs/$a/mat_$a.fa $dir/Chrs/$a/pat_$a.fa $snp_mode $dir/Lastz/$a/$a.lastz.txt $a " > $a.svmu.sh
	# sh  $a.svmu.sh
	cd ..
done
~~~~~


### run Assemblytics to detect CNVs

~~~~~bash
data="Alignment_nomask/Mummer"
assemblytics="software/Genome/Assemblytics-1.2/scripts/Assemblytics"
cat ../chr.list|while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	[ -e $a.delta ] || ln -s $data/$a/$a.delta
	# Assemblytics delta output_prefix unique_length_required min_size max_size
	$assemblytics $a.delta $a 10000 50 100000 && echo "$a done!"
	software/Genome/Assemblytics-1.2/scripts/Assemblytics_dotplot.R $a
	python3 software/Genome/Assemblytics-1.2/scripts/Assemblytics_within_alignment.py --delta $a.delta --min 50 --output $a.variants_within_alignments.bed
	cat $a.variants_within_alignments.bed $a.variants_between_alignments.bed|sort -k2n > $a.Assemblytics_structural_variants.bed

	cd ..
done

# filter by length and N rate
ref="Alignment_nomask/Chrs/$a"
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

cat ../list |while read a;do cd $a; python3 ../calCNVdensity.py $a.cnv.fil.bed chr.curated.len 1000000 1000000 $a >$a.cnv.1m.density;cd ..;done

# note: chr.curated.len is table format file that contains chromosome length.
~~~~~




### verify SVs

To generate a high-quality SV dataset, we manually checked all inversions and translocations with the following steps: a) clip 300 bp of upstream/downstream flanking sequence of each breakpoint between the two haplotypes, blast against local Pacbio reads with threshold identity >96% and aligned length >550 bp, and require the SV region where the maternal and paternal sequences aligned to have high similarity (>90%); b) if (a) failed, then check the 10X linked-read count between a 5 kb flanking region; c) if any breakpoint is not supported by 10X linked-reads, check the Hi-C heatmap of this region; if it shows an inversion or translocation pattern on heatmap or an ambiguous situation, then remove it.





<style class="fallback">body{visibility:hidden}</style><script>markdeepOptions={tocStyle:'long'};</script>
<!-- Markdeep: --><script src="https://casual-effects.com/markdeep/latest/markdeep.min.js?" charset="utf-8"></script>