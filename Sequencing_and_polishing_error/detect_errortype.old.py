#!/usr/bin/env python3

import sys
if len(sys.argv) < 2:
    print("Usage: python3 {} {} ".format(sys.argv[0], "*.base.xls.check"))
    exit()

def get_info(arr):
    """
    #CHROME    R_POS   R_REF   Q_POS   Q_REF   VAR_MUMMER  VAR_GATK    VAR_BCF NOTE_RG NOTE_DEP
    Chr Pos Depth   A   T   G   C   Pi  Chr Pos Depth   A   T   G   C   Pi  Chr Pos Depth   A   T
    G   C   Pi
    """
    matref = arr[2]
    patref = arr[4]
    tenx_bases = arr[13:18]
    mat_bases = arr[21:26]
    pat_bases = arr[-5:]

    return (matref, patref, tenx_bases, mat_bases, pat_bases)

def call_haplotype(base_dis, ty=None):
    bases = {'A': int(base_dis[0]),
             'T': int(base_dis[1]),
             'G': int(base_dis[2]),
             'C': int(base_dis[3])}
    pi = float(base_dis[4])
    sorted_bases = sorted(bases, key=lambda k: bases[k], reverse=True)
    hap1 = sorted_bases[0]
    hap2 = sorted_bases[1]
    if ty:
        if bases[hap2] == 0:
            hap2 = hap1 #this is homozygous
            return hap1, hap2
        elif pi < 0.25 and bases[hap2] >= 3: # the second most high depht base is not right;
            return False
        else:
            return hap1, hap2
    else:
        return hap1, hap2




with open(sys.argv[1], 'r') as fh:
    code = {'A': 2, 'T': 3, 'G': 5, 'C': 7}
    for i in fh:
        i = i.strip()
        arr = i.split("\t")
        if arr[20] == "-" or arr[-6] == "-":
            print(i + "\t" + "no_pacbio_read")
            continue
        (matref, patref, tenx_bases, mat_bases, pat_bases) = get_info(arr)
        assembly_pattern = code[matref] * code[patref]
        mat_hap_cor, mp = call_haplotype(mat_bases)
        pat_hap_cor, pp = call_haplotype(pat_bases)
        if call_haplotype(tenx_bases, 'tenx') == False:
            print(i + "\t" + "hard_to_identify")
            continue
        else:
            tenx_hap1, tenx_hap2 = call_haplotype(tenx_bases)
        tenx_pattern = code[tenx_hap1] * code[tenx_hap2]
        #print("tenx: {} {}".format(tenx_hap1, tenx_hap2))
        #print("assembly: {} {}".format(matref, patref))

        if assembly_pattern == tenx_pattern:
            if arr[6] == "-" and matref == mat_hap_cor and patref == pat_hap_cor:
                print(i + "\t" + "gatk_false_negative")
            else:
                print(i + "\t" + "right_polished")
        elif matref == tenx_hap1 and patref != tenx_hap2 and patref == pat_hap_cor:
            print(i + "\t" + "pat_sequencing_err")
        elif matref == tenx_hap1 and patref != tenx_hap2 and patref != pat_hap_cor:
            print(i + "\t" + "pat_polishing_err")
        elif matref == tenx_hap2 and patref != tenx_hap1 and patref == pat_hap_cor:
            print(i + "\t" + "pat_sequencing_err")
        elif matref == tenx_hap2 and patref != tenx_hap1 and patref != pat_hap_cor:
            print(i + "\t" + "pat_polishing_err")

        elif matref != tenx_hap1 and patref == tenx_hap2 and matref == mat_hap_cor:
            print(i + "\t" + "mat_sequencing_err")
        elif matref != tenx_hap1 and patref == tenx_hap2 and matref != mat_hap_cor:
            print(i + "\t" + "mat_polishing_err")
        elif matref != tenx_hap2 and patref == tenx_hap1 and matref == mat_hap_cor:
            print(i + "\t" + "mat_sequencing_err")
        elif matref != tenx_hap2 and patref == tenx_hap1 and matref != mat_hap_cor:
            print(i + "\t" + "mat_polishing_err")
        elif (matref != tenx_hap1 and patref != tenx_hap2) and (matref != tenx_hap2 and patref !=tenx_hap1):
            print(i + "\t" + "very_rare")
        else:
            print(i + "\t" + "other")

