#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

bedtools bamtobed -i $2 | \
mergeBed | \
bedtools complement -i - -g <(awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $1) > ${1%.*}.bed

samtools fasta $2 > ${2%.*}.fasta

cat <(awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $1) ${2%.*}.fasta > complete_set.fasta

samtools index $2

python tiling.py $1 $2 ${1%.*}.bed complete_set.fasta

printf "\nFinal tiles:\n"

cat ${1%.*}_decollapsed.agp

seqkit -is replace -p "^n+|n+$" -r "" ${1%.*}_decollapsed.fasta > tmp.fasta

cat tmp.fasta > ${1%.*}_decollapsed.fasta

rm tmp.fasta 