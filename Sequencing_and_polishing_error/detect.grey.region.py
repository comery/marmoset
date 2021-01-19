#!/usr/bin/env python3

import sys
if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "*.1coord", "pos-like.list"))
    exit()

aligned = {}
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split("\t")
        start = int(tmp[0])
        end = int(tmp[1])
        for b in range(start, end):
            aligned[b] = 1

with open(sys.argv[2], 'r') as fh:
    for i in fh:
        tmp = i.strip().split("\t")
        pos = tmp[1]
        if int(pos) not in aligned.keys():
            print(i.strip())


