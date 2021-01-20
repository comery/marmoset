#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    sys.exit('python3 %s <*.vcf.genotype>' % (sys.argv[0]))

inFile = sys.argv[1]

with open(inFile) as f:
    for line in f:
        line = line.rstrip()
        if "CHROM" in line:
            continue
        tmp = line.split('\t')
        temp = tmp[7].split(",")
        a1 = int(temp[0])
        a2 = int(temp[1])
        total = int(tmp[8])
        if a1 >0 and a2 >0 and total > 5 and total < 60:
            print(line)

