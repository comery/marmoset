#!/usr/bin/env python3

import sys
import pysam
from collections import defaultdict

def store(inFile):
	dic = defaultdict(dict)
	with open(inFile) as f:
		for line in f:
			if line[0] == '#': continue
			line = line.rstrip()
			tmp = line.split('\t')
			scaf = tmp[0]
			pos = int(tmp[1])
			ref = tmp[3]
			alt = tmp[4]
			dic[scaf][pos] = [ref, alt]
	return dic

def check(dic, bamFile):
	samfile = pysam.AlignmentFile(bamFile, "rb")

	for scaf in dic:
		for pos in sorted(dic[scaf].keys()):
			ref = dic[scaf][pos][0]
			alt = dic[scaf][pos][1]
			for pileupcolumn in samfile.pileup(scaf, pos-1, pos):
				if pileupcolumn.pos == pos-1:
					for pileupread in pileupcolumn.pileups:
						if not pileupread.is_del and not pileupread.is_refskip:
							flag = None
							if pileupread.alignment.query_sequence[pileupread.query_position] == ref:
								flag = 'ref'
							elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
								flag = 'alt'
							else:
								continue
							print('%s\t%i\t%s\t%s\t%s\t%s' % (scaf, pos, ref, alt, flag, pileupread.alignment.query_name))
	samfile.close()

def main():
	if len(sys.argv) != 3:
		sys.exit('python3 %s <tab> <bam>' % (sys.argv[0]))

	inFile = sys.argv[1]
	bamFile = sys.argv[2]

	dic = store(inFile)
	check(dic, bamFile)

if __name__ == '__main__':
	main()

