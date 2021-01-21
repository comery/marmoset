#!/usr/bin/env python3

import sys
sys.path.append('/hwfssz1/ST_DIVERSITY/PUB/USER/zhouyang/python_bin/util/')
from read_file import read_file
from collections import defaultdict
import numpy as np
from scipy.stats import ranksums

if len(sys.argv) != 5:
	sys.exit('python3 %s <bed> <matrix.add.filt> <sexLinked_HiC_scaffold.bed> <autosome.bed>' % (sys.argv[0]))

bedFile = sys.argv[1]
inFile = sys.argv[2]
sbedFile = sys.argv[3]
abedFile = sys.argv[4]

dic = {}
sdic = {}
adic = {}
interaction = defaultdict(lambda: defaultdict(list))

with open(bedFile) as f:
	for line in f:
		line = line.rstrip()
		tmp = line.split('\t')
		scaf = tmp[0]
		dic[scaf] = 1

with open(abedFile) as f:
	for line in f:
		line = line.rstrip()
		tmp = line.split('\t')
		scaf = tmp[0]
		adic[scaf] = 1

with open(sbedFile) as f:
	for line in f:
		line = line.rstrip()
		tmp = line.split('\t')
		scaf = tmp[0]
		sdic[scaf] = 1

with read_file(inFile) as f:
	for line in f:
		line = line.rstrip()
		tmp = line.split('\t')
		scaf0 = tmp[0]
		scaf1 = tmp[3]
		if tmp[6] == 'nan': continue
		strength = float(tmp[6])
		if scaf0 in dic:
			interaction[scaf0][scaf1].append(strength)
		if scaf1 in dic:
			interaction[scaf1][scaf0].append(strength)

for scaf in sorted(interaction.keys()):
	tmp_dic = {}
	n = 0
	q1 = 'NA'
	med = 'NA'
	q3 = 'NA'
	for scafID in sdic:
		if scafID in interaction[scaf]:
			if len(interaction[scaf][scafID]) > 5:
				arr = np.array(interaction[scaf][scafID])
				q1 = np.quantile(arr, 0.25)
				med = np.quantile(arr, 0.5)
				q3 = np.quantile(arr, 0.75)
				tmp_dic[scafID] = [arr, q1, med, q3]
	for k, v in sorted(tmp_dic.items(), key=lambda item: item[1][2], reverse=True):
		scafID = k
		arr = v[0]
		n = len(arr)
		q1 = v[1]
		med = v[2]
		q3 = v[3]
		break

	arr_dic = {}
	for scaf0 in interaction[scaf]:
		if scaf0 != scafID and scaf0 != scaf and scaf0 in adic:
			n0 = len(interaction[scaf][scaf0])
			#print(n0)
			#print(scaf0, scaf)
			#print(interaction[scaf][scaf0])
			if n0 >= 5:
				arr0 = np.array(interaction[scaf][scaf0])
				q10 = np.quantile(arr0, 0.25)
				med0 = np.quantile(arr0, 0.5)
				q30 = np.quantile(arr0, 0.75)
				arr_dic[scaf0] = [arr0, q10, med0, q30]
	if len(arr_dic) > 0:
		for k, v in sorted(arr_dic.items(), key=lambda item: item[1][2], reverse=True):
			scaf0 = k
			arr0 = v[0]
			q10 = v[1]
			med0 = v[2]
			q30 = v[3]
			n0 = len(arr0)
			if n > 0:
				s, pval = ranksums(arr, arr0)
			else:
				pval = 'NA'
	else:
		scaf0 = 'NA'
		n0 = 0
		q10 = 'NA'
		med0 = 'NA'
		q30 = 'NA'
		pval = 'NA'
	if q1 != 'NA': q1 = '%.4f' % (q1)
	if med != 'NA': med = '%.4f' % (med)
	if q3 != 'NA': q3 = '%.4f' % (q3)
	if q10 != 'NA': q10 = '%.4f' % (q10)
	if med0 != 'NA': med0 = '%.4f' % (med0)
	if q30 != 'NA': q30 = '%.4f' % (q30)
	if pval != 'NA': pval = '%.4f' % (pval)
	out = '\t'.join([scaf, scafID, str(n), q1, med, q3, scaf0, str(n0), q10, med0, q30, pval])
	print(out)

