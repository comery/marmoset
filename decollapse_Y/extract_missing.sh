#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

scaffold=$1

if [[ -e $scaffold/$scaffold.collapse/sda.asm.bam ]]; then

bedtools getfasta -fi $scaffold/$scaffold.fasta -bed <(	\
bedtools bamtobed -i $scaffold/$scaffold.collapse/sda.asm.bam | \
mergeBed | \
bedtools complement -i - -g <(awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $scaffold/$scaffold.fasta) | \
awk '{if ($3-$2>5000) print}') | \
seqkit -is replace -p "^n+|n+$" -r "" - > $scaffold/$scaffold.missing.fasta

if [[ $(seqtk comp $scaffold/$scaffold.missing.fasta | awk -v cutoff=$2 '{if ($9/$2>cutoff && $2-$9<20000) print $1}' | wc -l) > 0 ]]; then

	seqkit grep -v -f <(seqtk comp $scaffold/$scaffold.missing.fasta | awk -v cutoff=$2 '{if ($9/$2>cutoff && $2-$9<20000) print $1}') $scaffold/$scaffold.missing.fasta > $scaffold/$scaffold.missing.purged.fasta

	printf "%s: purged %s sequences.\n" $scaffold $(seqtk comp $scaffold/$scaffold.missing.fasta | awk -v cutoff=$2 '{if ($9/$2>cutoff && $2-$9<20000) print $1}' | wc -l)

	seqtk comp $scaffold/$scaffold.missing.fasta | awk -v cutoff=$2 '{if ($9/$2>cutoff&& $2-$9<20000) printf $0"\t"$9/$2"\n"}'
	
	mv $scaffold/$scaffold.missing.purged.fasta $scaffold/$scaffold.missing.fasta
else

	printf "%s: nothing to purge.\n" $scaffold

fi

else

	printf "%s: nothing to decollapse.\n" $scaffold

fi
