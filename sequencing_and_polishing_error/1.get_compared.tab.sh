cat list |while read a;do cd $a;echo "perl ../../bin/compare.snp.diff.pl $a.diff.snp $a.snp.vcf.gz $a.sambcf.vcf.gz $a $a.pileup $a.diff.1coords" >$a.work.sh; sel-qsub evo 20g 1 $a.work.sh;cd ..;done
