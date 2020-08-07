#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "*.primers.txt", "*.blast.good"))
    exit()

def get_info(record):
    tmp = record.strip().split("\t")
    info = tmp[0].split("_")
    chrome = info[0]
    target_pos = info[-2]
    art_id = chrome + "_" + target_pos
    return art_id


def main():
    pos_len = {}
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split("\t")
            a = tmp[0] + "_" + tmp[1]
            pos_len[a] = int(tmp[2])
    c1 = open("class1.primer.txt", 'w')
    c2 = open("class2.primer.txt", 'w')
    c3 = open("class3.primer.txt", 'w')
    c4 = open("class4.primer.txt", 'w')
    c5 = open("class5.primer.txt", 'w')
    c6 = open("class6.primer.txt", 'w')

    with open(sys.argv[1], 'r') as fh:
        good_for = set()
        good_rev = set()
        good_both = set()
        for i in fh:
            if i.startswith("Index"):
                continue
            art_id = get_info(i)
            if pos_len[art_id] == 1:
                c1.write(i)
            elif pos_len[art_id] > 1 and pos_len[art_id] <= 10:
                c2.write(i)
            elif pos_len[art_id] > 10 and pos_len[art_id] <= 20:
                c3.write(i)
            elif pos_len[art_id] > 20 and pos_len[art_id] <= 30:
                c4.write(i)
            elif pos_len[art_id] > 30 and pos_len[art_id] <= 40:
                c5.write(i)
            elif pos_len[art_id] > 40 and pos_len[art_id] <= 50:
                c6.write(i)
    c1.close()
    c2.close()
    c3.close()
    c4.close()
    c5.close()
    c6.close()


if __name__ == '__main__':
    main()
