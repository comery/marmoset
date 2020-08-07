# get a potential list of ins
#awk '$2==$3 && $8-$7>200 && $8-$7<500' Chr1.ins.tab |head -50 >ins.potential.list
#awk '$3-$2>200 && $3-$2<500 && $7==$8' Chr1.del.tab |head -50 >del.potential.list
python3 extractTemplate.py ins.potential.list mat_Chr1.fa pat_Chr1.fa
python3 extractTemplate.py del.potential.list mat_Chr1.fa pat_Chr1.fa
awk '$2==$3 && $8-$7>200 && $8-$7<500' Chr1.ins.tab |sed -n '51,$p' >ins.potential2.list
awk '$3-$2>200 && $3-$2<500 && $7==$8' Chr1.del.tab |sed -n '51,$p' >del.potential2.list
python3 extractTemplate.py ins.potential2.list mat_Chr1.fa pat_Chr1.fa
python3 extractTemplate.py del.potential2.list mat_Chr1.fa pat_Chr1.fa
