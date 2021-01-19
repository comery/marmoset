#!/usr/bin/env python3
import os
import sys
import subprocess
if len(sys.argv) < 2:
    print("Usage: python3 {} {}".format(sys.argv[0], "sv.tab"))
    exit()


def main():

    with open(sys.argv[1], 'r') as fh, open("err.log", 'w') as log:
        for r in fh:
            if r.startswith("#"):
                continue
            tmp = r.strip().split("\t")
            svtype = tmp[-2]
            mat_name = tmp[0]
            mat_start = int(tmp[1])
            mat_end = int(tmp[2])
            pat_name = tmp[5]
            pat_start = int(tmp[6])
            pat_end = int(tmp[7])

            # all positions
            if mat_start == mat_end:
                if os.path.exists("INS") == False:
                    os.mkdir("INS")
                sv_len = pat_end - pat_start + 1
                # extend to 1000 bp 
                shift = int(500 - (sv_len)/2)
                mat_extend_f = mat_start - shift
                mat_extend_r = mat_end + shift
                pat_extend_f = pat_start - shift
                pat_extend_r = pat_end + shift
                print("{}_{}_{}\t{}".format(pat_name, pat_extend_f, pat_extend_r, sv_len))
            else:
                if os.path.exists("DEL") == False:
                    os.mkdir("DEL")
                sv_len = mat_end - mat_start + 1
                shift = int(500 - (sv_len)/2)
                mat_extend_f = mat_start - shift
                mat_extend_r = mat_end + shift
                pat_extend_f = pat_start - shift
                pat_extend_r = pat_end + shift
                print("{}_{}_{}\t{}".format(mat_name, mat_extend_f, mat_extend_r, sv_len))


if __name__ == '__main__':
    main()
