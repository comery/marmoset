#!/bin/bash

#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

scaffold=$1

if [[ -s $scaffold/$scaffold.collapse/sda.asm.bam && $(samtools view -c $scaffold/$scaffold.collapse/sda.asm.bam) > 0 ]]; then

cd $scaffold

mkdir -p tiling

cd tiling

ln -sf ../$scaffold.fasta
ln -sf ../$scaffold.collapse/sda.asm.bam

ln -sf ../../../tiling/tiling.sh
ln -sf ../../../tiling/tiling.py

bash tiling.sh $scaffold.fasta sda.asm.bam

cat ${scaffold}_decollapsed.fasta >> ../../chrY_unplaced_decollapsed.fasta
cat ${scaffold}_decollapsed.agp >> ../../chrY_unplaced_decollapsed.agp
cat ${scaffold}_unplaced.fasta >> ../../chrY_unplaced_SDA_unplaced.fasta

cd ../..

else

	printf "%s: no decollapsed sequence.\n" $scaffold

fi
