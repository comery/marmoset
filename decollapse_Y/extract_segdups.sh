#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

scaffold=$1

if [[ -e $scaffold/$scaffold.collapse/sda.asm.bam ]]; then

samtools fasta $scaffold/$scaffold.collapse/sda.asm.bam 2> /dev/null | awk 'NR==FNR{a[$1]=$2;next} {for(k in a) gsub(k,a[k]);print}' <(bedtools bamtobed -i $scaffold/$scaffold.collapse/sda.asm.bam | awk '{printf $4"\t"$1":"$2"-"$3"_"$4"\n"}') - > $scaffold/$scaffold.segdups.fasta

else

	printf "%s: no decollapsed sequence.\n" $scaffold

fi
