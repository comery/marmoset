#!/usr/bin/env python3

import sys

if len(sys.argv) < 2:
    sys.exit("Usage: python3 {} {}".format(sys.argv[0], "*.primers"))

def main():
    short = {"FORWARD": 'F',
             'REVERSE': 'R'}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            i= i.strip()
            if i.startswith("Index"):
                continue
            tmp = i.split("\t")
            info = tmp[1]
            index = tmp[2]
            strand = tmp[4]
            start = int(tmp[5])
            prilen = int(tmp[6])
            pri = tmp[11]
            tmp1 = info.split("_")
            target = tmp1[-1]
            chrom = tmp1[0] + "_" + tmp1[1]
            template_start = int(tmp1[2].split("-")[0])
            pri_start = template_start + (start -1)
            if short[strand] == "F":
                pri_end = pri_start + (prilen - 1)
            else:
                pri_end = pri_start - (prilen - 1)

            print(">{}_{}_{}_{}_{}_{}\n{}".format(chrom, target, index, short[strand], pri_start, pri_end,
                                                 pri))

if __name__ == '__main__':
    main()
