#!/usr/bin/env python3

import sys
import pysam
from collections import defaultdict

def examine(inFile, bamFile):
	dic = defaultdict(dict)
	previous = ''
	current = ''

	with open(inFile) as f:
		for line in f:
			if line[0] == '#': continue
			line = line.rstrip()
			if previous == '':
				previous = line
				continue
			else:
				current = line
				
			tmp1 = previous.split('\t')
			scaf1 = tmp1[0]
			pos1 = int(tmp1[1])
			(ref1, alt1) = tmp1[6].split('|')[:]

			tmp2 = current.split('\t')
			scaf2 = tmp2[0]
			pos2 = int(tmp2[1])
			(ref2, alt2) = tmp2[6].split('|')[:]
			
			if scaf1 == scaf2:	
				combination = check(previous, current, bamFile)

				out = '%s\t%i\t%s/%s\t%i\t%s/%s\t' % (scaf1, pos1, ref1, alt1, pos2, ref2, alt2)
				#print(previous)
				#print(current)
				for i in sorted(combination.keys()):
					out += '%s\t%i\t' % (i, combination[i])
				out = out.rstrip()
				print(out)
			previous = line
	
	return dic

def extract_reads(bamFile, scaf, pos, ref, alt):
	samfile = pysam.AlignmentFile(bamFile, 'rb')
	reads = {}
	for pileupcolumn in samfile.pileup(scaf, pos-1, pos):
		if pileupcolumn.pos == pos-1:
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					flag = None
					if pileupread.alignment.query_sequence[pileupread.query_position] == ref:
						reads[pileupread.alignment.query_name] = 'ref'
					elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
						reads[pileupread.alignment.query_name] = 'alt'
					else:
						continue
	samfile.close()
	return reads

def check(previous, current, bamFile):
	tmp1 = previous.split('\t')
	scaf1 = tmp1[0]
	pos1 = int(tmp1[1])
	(ref1, alt1) = tmp1[6].split('|')[:]

	tmp2 = current.split('\t')
	scaf2 = tmp2[0]
	pos2 = int(tmp2[1])
	(ref2, alt2) = tmp2[6].split('|')[:]

	if scaf1 == scaf2:
		reads1 = extract_reads(bamFile, scaf1, pos1, ref1, alt1)
		reads2 = extract_reads(bamFile, scaf2, pos2, ref2, alt2)
		
		combination = {'ref-ref': 0, 'ref-alt': 0, 'alt-ref': 0, 'alt-alt': 0}
		for id in reads1:
			if id in reads2:
				combination[reads1[id] + '-' + reads2[id]] += 1
		return combination

def main():
	if len(sys.argv) != 3:
		sys.exit('python3 %s <genotype_genomicDBI_gather.g.vcf.tab.snp.mat.filt> <mother.bam>' % (sys.argv[0]))

	inFile = sys.argv[1]
	bamFile = sys.argv[2]

	examine(inFile, bamFile)

if __name__ == '__main__':
	main()

