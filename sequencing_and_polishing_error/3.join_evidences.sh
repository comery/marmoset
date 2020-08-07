cat list |while read a
do 
	cd $a
	#python3 ../../bin/join_5evidence.py $a.compared.tab $a.rqpos.coords $a.10x.compared.base.xls $a.mat2mat.pb.cor.compared.base.xls $a.mat2mat.pb.raw.compared.base.xls $a.pat2pat.pb.cor.compared.base.xls $a.pat2pat.pb.raw.compared.base.xls >$a.5evidence.xls
	python3 ../../bin/snp_5evidence.statistic.v2.py $a.5evidence.xls
	grep "True" $a.5evidence.xls.tag|cut -f 1-8 > SNP.checked.true.xls
	cd ..
done

# all TRUE SNP (10x confirmed)
cat list|while read a;do cut -f 1-5,7 $a/SNP.checked.true.xls ;done > all.SNP.checked.true.xls

# error type statistic
cat list|while read a;do perl ../bin/error.stat.pl $a/$a.5evidence.xls.tag ;done >all.errortype.summary.xls

# the intersection of GATK and Samtools
cat list|while read a; do grep -v '^#' $a/$a.5evidence.xls|cut -f 5|grep -v '-'|wc -l; done
