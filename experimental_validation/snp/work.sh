# generate template fasta base SNP positions
a1.extracFlankingSeq.sh
primerDesign.snp.template.fa

# submit these template fasta upto web tools
http://batchprimer3.bioinformatics.ucdavis.edu/cgi-bin/batchprimer3/batchprimer3.cgi?PRIMER_TYPE=4&CLEAR_FORM=no&LOAD_EXAMPLE=no

# download primer results and collect into all
batchPrimer3.primer.txt

# convert primers.txt to primer.fa
python3 ../changePrimerTab2Fa.py batchPrimer3.primer.txt > batchPrimer3.primer.fa

# conduct primer-blast using blastn
blastn -query batchPrimer3.primer.fa -out batchPrimer3.primer.blastout -db /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/SNP_INDEL/experiments/blastdb/maternal -outfmt 6  -max_target_seqs 5 -num_threads 8 -task blastn

# find specific primers, (identity <80%), output primer ids
python3 ../check_PrimerBlastResult.py batchPrimer3.primer.blastout 0.8 >batchPrimer3.primer.blastout.good
sort batchPrimer3.primer.blastout.good > batchPrimer3.primer.blastout.good.sort
# retrieve primers based on primer id
python3 ../pickGoodPrimers.py batchPrimer3.primer.txt batchPrimer3.primer.blastout.good.sort

# format primers
python3 ../format_primer.py batchPrimer3.primer.blastout.good.sort.primers 1 >both_specific.good_top1.primer.txt
