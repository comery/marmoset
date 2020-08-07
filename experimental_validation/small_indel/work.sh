blastn -query all.primer.fa -out all.primer.blastout -db /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/SNP_INDEL/experiments/blastdb/maternal -outfmt 6  -max_target_seqs 5 -num_threads 8 -task blastn


awk '($6>100 && $6<200) || ($6>350 && $6<450)' all.primer.blastout.good.sort.primers
