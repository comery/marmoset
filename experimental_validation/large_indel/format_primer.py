#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "*.primers.txt", "index threshold"))
    exit()


def parser_blast(file):
    pre = ""
    target_match = []
    with open(file, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if not line:
                break
            tmp = line.strip().split("\t")
            if pre == "":
                target_match = [line,]
                pre = tmp[1]
            elif tmp[1] == pre:
                target_match.append(line)
            else:
                yield target_match
                target_match = [line,]
                pre = tmp[1]
        yield target_match


def main():
    primers = {}
    unique = {}
    cutoff = int(sys.argv[2])
    for i in parser_blast(sys.argv[1]):
        for r in i:
            tmp = r.split("\t")
            index = tmp[2]
            if int(index) > 5:
                continue
            if int(index) > cutoff:
                continue
            pri_id = tmp[1] + "_" + index
            strand = tmp[4]
            pri = tmp[11]
            if tmp[1] in unique.keys():
                continue
            else:
                if pri_id in primers.keys():
                    primers[pri_id].append(strand + "\t" + pri)
                    if len(primers[pri_id]) == 2:
                        print("\t".join(primers[pri_id]))
                        unique[tmp[1]] = 1
                else:
                    primers[pri_id] = [pri_id + "\t" + strand + "\t" + pri, ]


if __name__ == '__main__':
    main()
