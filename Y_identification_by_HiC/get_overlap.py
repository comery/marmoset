#!/usr/bin/env python3

import sys
sys.path.append('/hwfssz1/ST_DIVERSITY/PUB/USER/zhouyang/python_bin/util/')
from read_file import read_file
from collections import defaultdict

if len(sys.argv) != 3:
	sys.exit("python %s <bed> <matrix.add>" % (sys.argv[0]))
script, bed1, bed2 = sys.argv


dic = defaultdict(list)
f0 = open(bed1, "r")
for line in f0:
	line = line.strip().split()
	ID = line[0]
	bg = int(line[1])
	ed = int(line[2])
	dic[ID].append([bg, ed])

with read_file(bed2) as f1:
	for line in f1:
		line = line.strip()
		tmp = line.split()
		scaf0 = tmp[0]
		bg0 = int(tmp[1])
		ed0 = int(tmp[2])
		scaf1 = tmp[3]
		bg1 = int(tmp[4])
		ed1 = int(tmp[5])
		if scaf0 in dic:
			for i in dic[scaf0]:
				bg, ed = i[:]
				if bg0 >= bg and ed0 <= ed:
					print(line)
		elif scaf1 in dic:
			for i in dic[scaf1]:
				bg, ed = i[:]
				if bg1 >= bg and ed1 <= ed:
					print(line)

