#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "*.primers", "*.blast.good"))
    exit()

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
            tmp = i.strip().split("\t")
            forward_id = tmp[0] + "_" + tmp[3] + "_" + tmp[4]
            reverse_id = tmp[0] + "_" + tmp[8] + "_" + tmp[9]
            if forward_id in good_primers and reverse_id not in good_primers:
                good_for.add(tmp[1])
            elif forward_id not in good_primers and reverse_id in good_primers:
                good_rev.add(tmp[1])
            elif forward_id in good_primers and reverse_id in good_primers:
                good_both.add(tmp[1])
                print(i.strip(), file=out)
    out.close()

    print("good_for:\t{}".format(len(good_for)))
    print("good_rev:\t{}".format(len(good_rev)))
    print("good_both:\t{}".format(len(good_both)))



if __name__ == '__main__':
    main()
