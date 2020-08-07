#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "<pileup file>", "<rqpos.coords>"))
    exit()

def get_pos(file):
    pos = {}
    with open(file ,'r') as fh:
        for i in fh.readlines():
            if i.startswith("#"):
                continue
            tmp = i.strip().split("\t", 3)[2]
            if tmp == "-" or tmp == ".":
                continue
            else:
                pos[int(tmp)] = 1
    return pos


def main():
    pos = get_pos(sys.argv[2])
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            line = line.strip()
            loca = line.split("\t", 2)[1]
            if int(loca) in pos.keys():
                print(line)


if __name__ == '__main__':
    main()

