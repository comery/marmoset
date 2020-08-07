#!/usr/bin/env python3
import os
import sys
import subprocess
if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {}".format(sys.argv[0], "sv.tab",
                                                   "mat.fa", "pat.fa"))
    exit()


def read_one_fa(fa):
    name = ""
    seq = ""
    with open(fa, 'r') as fh:
        for i in fh:
            if i.startswith(">"):
                name = i.strip().replace(">", "")
            else:
                seq += i.strip()
    return name, seq

def out2fa(f, string):
    with open(f, 'w') as fw:
        fw.write(string)


def reverse_complementation(sequence):
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence[::-1]


def aln_shell(outdir):
    with open(outdir + "/work.sh", 'w') as fw:
        fw.write("cat mat.fa pat.fa >mat.pat.fa\n")
        fw.write("mafft --thread 2 --quiet mat.pat.fa >mat.pat.aln \n")


def main():
    mat_name, matref = read_one_fa(sys.argv[2])
    pat_name, patref = read_one_fa(sys.argv[3])
    matref_len = len(matref)
    patref_len = len(patref)

    with open(sys.argv[1], 'r') as fh, open("err.log", 'w') as log:
        extend = 300
        for r in fh:
            if r.startswith("#"):
                continue
            tmp = r.strip().split("\t")
            svtype = tmp[-2]
            mat_start = int(tmp[1])
            mat_end = int(tmp[2])
            pat_start = int(tmp[6])
            pat_end = int(tmp[7])

            # all positions
            if mat_start == mat_end:
                if os.path.exists("INS") == False:
                    os.mkdir("INS")
                outdir = "INS/{}_{}".format(mat_start, mat_end)
                sv_len = pat_end - pat_start + 1
                # extend to 1000 bp 
                shift = int(500 - (sv_len)/2)
                mat_extend_f = mat_start - shift
                mat_extend_r = mat_end + shift
                pat_extend_f = pat_start - shift
                pat_extend_r = pat_end + shift
            else:
                if os.path.exists("DEL") == False:
                    os.mkdir("DEL")
                outdir = "DEL/{}_{}".format(pat_start, pat_end)
                sv_len = mat_end - mat_start + 1
                shift = int(500 - (sv_len)/2)
                mat_extend_f = mat_start - shift
                mat_extend_r = mat_end + shift
                pat_extend_f = pat_start - shift
                pat_extend_r = pat_end + shift
            if os.path.exists(outdir) == False:
                os.mkdir(outdir)
            mat_seq = matref[mat_extend_f-1:mat_extend_r]
            pat_seq = patref[pat_extend_f-1:pat_extend_r]
            out2fa(outdir + "/mat.fa", ">{}_{}_{}\n{}\n".format(mat_name,
                                                               mat_extend_f,
                                                               mat_extend_r,
                                                               mat_seq))
            out2fa(outdir + "/pat.fa", ">{}_{}_{}\n{}\n".format(pat_name,
                                                                pat_extend_f,
                                                                pat_extend_r,
                                                                pat_seq))

            aln_shell(outdir)
            wkdir = os.getcwd()
            os.chdir(outdir)
            subprocess.call("sh work.sh 1>>std.log 2>>std.err", shell=True)
            os.chdir(wkdir)


if __name__ == '__main__':
    main()
