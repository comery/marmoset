#!/usr/bin/env python3
import os
import sys
import subprocess
if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {} ".format(sys.argv[0], "Assemblytics_structural_variants.bed",
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

def gc_content(sequence):
    sequence = sequence.upper()
    seqlen = len(sequence)
    a = sequence.count("A")
    t = sequence.count("T")
    c = sequence.count("C")
    g = sequence.count("G")
    n = sequence.count("N")
    gc_content = (g + c) / seqlen
    n_content = n / seqlen
    return gc_content, n_content

def get_pos(info):
    tmp = info.split(":")
    tmp1 = tmp[1].split("-")
    return tmp[0], int(tmp1[0]), int(tmp1[1])

def main():
    mat_name, matref = read_one_fa(sys.argv[2])
    pat_name, patref = read_one_fa(sys.argv[3])
    matref_len = len(matref)
    patref_len = len(patref)
    with open(sys.argv[1], 'r') as fh:
        for r in fh:
            if r.startswith("#"):
                continue
            elif 'Insertion' in r or 'Deletion' in r:
                continue
            tmp = r.strip().split("\t")
            svtype = tmp[6]
            ref_id = tmp[0]
            mat_start = int(tmp[1])
            mat_end = int(tmp[2])
            mat_varlen = abs(mat_end - mat_start)
            que_id, pat_start, pat_end = get_pos(tmp[9])
            pat_varlen = abs(pat_end - pat_start)
            varseq_mat = matref[mat_start-1:mat_end]
            varseq_pat = patref[pat_start-1:pat_end]
            mat_gc_content, mat_n_content = gc_content(varseq_mat)
            pat_gc_content, pat_n_content = gc_content(varseq_pat)

            var_rate = "{:.4f}".format(pat_varlen/mat_varlen)

            print("{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{}".format(ref_id, mat_start, mat_end,
                                                                                              mat_varlen,
                                                                                      mat_gc_content, mat_n_content,
                                                                                      que_id, pat_start, pat_end,
                                                                                              pat_varlen,
                                                                                      pat_gc_content, pat_n_content,
                                                                                      var_rate))


if __name__ == '__main__':
    main()
