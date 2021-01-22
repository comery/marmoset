#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

for scaffold in *.1; do 
if [ -f "$scaffold/$scaffold.collapse/sda.assemblies.fasta" ]; then 
printf "%s\t%s\t%s\n" 	$scaffold \
			$(grep -c ">" $scaffold/$scaffold.collapse/sda.assemblies.fasta) \
			$(samtools view -c $scaffold/$scaffold.collapse/sda.asm.bam)

else

printf "%s: no decollapse occurred.\n"   $scaffold

fi

done
