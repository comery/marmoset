#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    print("Usage: python3 {} {}".format(sys.argv[0], "*.blast"))
    exit()

def main():
    bad_list = set()
    all_list = set()
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.rstrip().split()
            all_list.add(tmp[0])
            tmp1 = tmp[0].split("_")
            query_chr = tmp1[0] + "_" + tmp1[1]
            target_chr = tmp[1]
            prilen = abs(int(tmp1[-1]) - int(tmp1[-2])) + 1
            # itself
            if query_chr == target_chr and tmp1[-2] == tmp[8] and tmp1[-1] == tmp[9]:
                continue
            else:
                if tmp[0] in bad_list:
                    continue
                aligned = int(tmp[3])
                mismatched = int(tmp[4])
                gap = int(tmp[5])
                matched = aligned - mismatched - gap
                real_id = matched / prilen
                three_end_tense = prilen - int(tmp[7])
                if real_id > 0.8:
                    bad_list.add(tmp[0])

    good_list = all_list - bad_list
    with open(sys.argv[1] + ".good", 'w') as fh:
        for i in good_list:
            print(i, file=fh)

if __name__ == '__main__':
    main()
