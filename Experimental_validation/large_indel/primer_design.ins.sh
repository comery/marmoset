# generate template fasta base SNP positions
cd ../
python3 extractTemplate.py ins.potential.list mat_Chr1.fa pat_Chr1.fa
cd INS
cat */pat.fa > large.ins.template.fa

# submit these template fasta upto web tools
http://batchprimer3.bioinformatics.ucdavis.edu/cgi-bin/batchprimer3/batchprimer3.cgi?PRIMER_TYPE=4&CLEAR_FORM=no&LOAD_EXAMPLE=no

# download primer results and collect into all
large.ins.primer.txt

# convert primers.txt to primer.fa
python3 ../changePrimerTab2Fa.py large.ins.primer.txt > large.ins.primer.fa

# conduct primer-blast using blastn
blastn -query large.ins.primer.fa -out large.ins.primer.blastout -db /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/SNP_INDEL/experiments/blastdb/paternal -outfmt 6  -max_target_seqs 5 -num_threads 4 -task blastn

# find specific primers, (identity <80%), output primer ids
python3 ../check_PrimerBlastResult.py large.ins.primer.blastout 0.8 >large.ins.primer.blastout.good
sort large.ins.primer.blastout.good > large.ins.primer.blastout.good.sort
# retrieve primers based on primer id
python3 ../pickGoodPrimers.py large.ins.primer.txt large.ins.primer.blastout.good
# it will get batchPrimer3.large.del.primer.blastout.good.sort.primers

# format primers
python3 ../format_primer.py large.ins.primer.blastout.good.sort.primers  100 > chr1_large_ins.primer.txt



