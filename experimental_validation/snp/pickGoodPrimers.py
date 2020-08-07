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
    art_id = chrome + "_" + target_pos + "_" + index + "_" + strand
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
            """
            if forward_id in good_primers and reverse_id not in good_primers:
                good_for.add(tmp[1])
            elif forward_id not in good_primers and reverse_id in good_primers:
                good_rev.add(tmp[1])
            elif forward_id in good_primers and reverse_id in good_primers:
                good_both.add(tmp[1])
                print(i.strip(), file=out)
            """
    out.close()

    print("good_for:\t{}".format(len(good_for)))
    print("good_rev:\t{}".format(len(good_rev)))
    print("good_both:\t{}".format(len(good_both)))



if __name__ == '__main__':
    main()
