#!/usr/bin/env python3

import sys
import argparse

inFile = sys.argv[1]

def chimeric(a1, a2):
    rate1 = 2 * a2 / (a1 + a2)
    rate2 = (a1 - a2) / (a1 + a2)
    rate3 = a2 / (a1 + a2)
    rate4 = (2 * a1 + a2) / (a2 -a1)
    rate5 = a1 / (a1 + a2)
    rate6 = 2 * a1 / (a1 + a2)

    a =  [rate1, rate2, rate3, rate4, rate5, rate6]
    b = ""
    for i in a:
        if i >0 and i < 0.5:
            b += str(i) + "\n"

    return b



with open(inFile) as f:
    for line in f:
        line = line.rstrip()
        tmp = line.split('\t')
        allele1 = tmp[2]
        allele2 = tmp[4]
        tmp1 = tmp[7].split(",")
        allele1_count = int(tmp1[0])
        allele2_count = int(tmp1[1])
        total_count = int(tmp[8])
        if allele1_count == 0 or allele2_count == 0:
            continue
        elif total_count < 5 or total_count > 60:
            continue
        elif allele1_count/allele2_count >=0.95 and allele1_count/allele2_count <1.05:
            continue
        elif allele2_count/allele1_count >=0.95 and allele2_count/allele1_count <1.05:
            continue
        #print(line + "\t" + chimeric(allele1_count, allele2_count))
        print( chimeric(allele1_count, allele2_count))
