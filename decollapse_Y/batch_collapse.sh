#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

while read -r line;
do

mkdir collapse/${line}

samtools faidx GCA_011100555.1_mCalJac1.pat.X_genomic.reorderedY.fasta ${line} > collapse/${line}/${line}.fasta
samtools faidx collapse/${line}/${line}.fasta

samtools view -b noAutosomes.noX.bam ${line} > collapse/${line}/${line}.bam
samtools index collapse/${line}/${line}.bam

sbatch -p vgl,hpc --cpus-per-task=24 --output=logs/${line}_%j.log --error=logs/${line}_%j.log SDA.sh collapse/${line}/${line}.fasta collapse/${line}/${line}.bam 33.5 collapse/${line}/${line}.collapse

done < ${1}
