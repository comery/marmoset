# get a potential list of ins
awk '$2==$3 && $8-$7>500 && $8-$7<1000' Chr1.ins.tab |head -50 >ins.potential3.list
awk '$3-$2>500 && $3-$2<1000 && $7==$8' Chr1.del.tab |head -50 >del.potential3.list
python3 extractTemplate.py ins.potential3.list mat_Chr1.fa pat_Chr1.fa
python3 extractTemplate.py del.potential3.list mat_Chr1.fa pat_Chr1.fa
