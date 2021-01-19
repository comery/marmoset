#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "*.blast", "threshold"))
    exit()

def main():
    threshold = float(sys.argv[2])
    bad_list = set()
    all_list = set()
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.rstrip().split()
            tmp1 = tmp[0].split("_")
            pri_sign = "_".join(tmp1[1:6])
            all_list.add(pri_sign)
            query_chr = tmp1[1]
            target_chr = tmp[1].split("_")[-1]
            prilen = abs(int(tmp1[-1]) - int(tmp1[-2])) + 1
            # itself
            shift_s = abs(int(tmp1[-2]) - int(tmp[8]))
            shift_e = abs(int(tmp1[-1]) - int(tmp[9]))
            if query_chr == target_chr and shift_s < 10 and shift_e < 10:
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
                if real_id > threshold:
                    bad_list.add(pri_sign)

    good_list = all_list - bad_list
    for i in good_list:
        print(i)

if __name__ == '__main__':
    main()
