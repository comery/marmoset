#!/usr/bin/env python3
import sys

if len(sys.argv) < 5:
    print("Usage: python3 {} {} {} {} {}".format(sys.argv[0],
                                                 "All.chr.diff.snp",
                                                 "chr.curated.len",
                                                 '<win length>',
                                                 '<step>',
                                                 "<Chr[int]>") )
    exit()

lens = {}

with open(sys.argv[2], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        lens[tmp[0]] = int(tmp[1])

def pos2dict(fh):
    pos = {}
    for i in fh:
        if i.startswith("#"):
            continue
        tmp = i.strip().split()
        if "-" in tmp[1]:
            mat_pos = int(tmp[1].split("-")[0])
        elif "_" in tmp[1]:
            mat_pos = int(tmp[1].replace("_", ""))
        else:
            mat_pos = int(tmp[1])

        pos[mat_pos] = 1
    return pos

def count_var(pos, start, end):
    vars = 0
    for i in range(start, end):
        if i in pos.keys():
            vars += 1

    return vars


chrom = "mat_" + sys.argv[5]

chrom_len = lens[chrom]

with open(sys.argv[1], 'r') as sp:
    start = 0
    pos = pos2dict(sp)
    win = int(sys.argv[3])
    step = int(sys.argv[4])
    end = start + win
    while end <= chrom_len:
        vars = count_var(pos, start, end)
        rate = vars / win
        print("{}\t{}\t{}\t{}\t{:.5f}".format(sys.argv[4], start, end, vars, rate))
        start += step
        end = start + win

    vars = count_var(pos, start, chrom_len)
    rate = vars/ (chrom_len - start + 1)
    print("{}\t{}\t{}\t{}\t{:.5f}".format(sys.argv[4], start, chrom_len, vars, rate))



