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
            if tmp[0] == tmp[1] + "_" + tmp[8] + "_" + tmp[9]:
                continue
            else:
                if tmp[0] in bad_list:
                    continue
                tmp1 = tmp[0].split("_")
                querry = tmp1[0] + "_" + tmp1[1]
                target = tmp[1]
                identity = float(tmp[2])
                if identity >= 80:
                    bad_list.add(tmp[0])

    good_list = all_list - bad_list
    with open(sys.argv[1] + ".good", 'w') as fh:
        for i in good_list:
            print(i, file=fh)

if __name__ == '__main__':
    main()
