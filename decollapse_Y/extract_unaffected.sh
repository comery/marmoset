#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

rm -f all_unaffected.fasta

for scaffold in *.1; do

if ! [[ -s $scaffold/$scaffold.collapse/sda.asm.bam && $(samtools view -c $scaffold/$scaffold.collapse/sda.asm.bam)>0 ]]; then

cat $scaffold/$scaffold.fasta >> all_unaffected.fasta

else

	printf "%s: sequence decollapsed. skipping.\n" $scaffold

fi

done
