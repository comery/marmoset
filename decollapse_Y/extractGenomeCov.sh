#++++     Credit: Giulio Formenti giulio.formenti@gmail.com     ++++

#bed from scaffold list
grep -f y.sda.ls mCalJac1.combine.sda.fasta.fai | awk '{print $1"\t1\t"$2}' > y.sda.bed

# view chr1 only
samtools view aligned_reads.bam -L y.sda.bed -@ 24 -b > chrYonly.sda.bam

#genomecov
bedtools genomecov -ibam chrYonly.sda.bam > chrYonly.sda.genomecov

#Total coverage
grep genome chrYonly.sda.genomecov

#Null bases to remove
grep -vf y.sda.ls chrYonly.sda.genomecov | grep -v genome | awk '{sum+=$3}END{print sum}'

#gaps to remove
grep -f  y.sda.ls mCalJac1.combine.sda.gaps | awk '{sum+=$2}END{print sum}'