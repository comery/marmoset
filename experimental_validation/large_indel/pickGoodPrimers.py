#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "*.primers.txt", "*.blast.good"))
    exit()

def get_info(record):
    short = {"FORWARD": 'F',
             "REVERSE": 'R'}
    tmp = record.strip().split("\t")
    info = tmp[1].split("_")
    chrome = info[0]
    target_pos = info[-1]
    index = tmp[2]
    strand = short[tmp[4]]
    art_id = "_".join(info[1:]) + "_" + index + "_" + strand
    return art_id


def main():
    good_primers = set()
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            good_primers.add(i.strip())
    out = open(sys.argv[2] + ".primers", 'w')
    with open(sys.argv[1], 'r') as fh:
        good_for = set()
        good_rev = set()
        good_both = set()
        for i in fh:
            if i.startswith("Index"):
                continue
            art_id = get_info(i)
            if art_id in good_primers:
                print(i.strip(), file=out)
    out.close()


if __name__ == '__main__':
    main()
