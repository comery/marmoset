#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    print("python3 {} snp.xls, maternal.allele.stat, paternal.allel.stat".format(sys.argv[0]))
    sys.exit()

def addtwodimdict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})


def minor_region(f):
    collapse = {}
    with open(f, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.rstrip().split("\t")
            if tmp[8] == "minor":
                snps = tmp[2]
                chrom = tmp[0]
                poses = snps.replace("(", "").replace(")", "").split(",")
                poses = list(map(lambda x:int(x), poses))
                for p in range(poses[0]-1,poses[-1]):
                    addtwodimdict(collapse, chrom, p, 1)
                    #collapse[p] = 1
    return collapse



def main():
    maternal_collapse = minor_region(sys.argv[2])
    paternal_collapse = minor_region(sys.argv[3])
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            i = i.strip()
            tmp = i.split("\t")
            chrom = tmp[0].replace("mat_", "")
            mat_pos1 = int(tmp[1])
            mat_pos2 = int(tmp[2])
            pat_pos1 = int(tmp[3])
            pat_pos2 = int(tmp[4])
            mat_sign = "-"
            mat_tag = 0
            pat_sign = "-"
            pat_tag = 0
            for m in range(mat_pos1, mat_pos2+1):
                if m in maternal_collapse[chrom].keys() and maternal_collapse[chrom][m]:
                    mat_tag += 1
            if mat_tag > 0:
                mat_sign = "minor"
            for p in range(pat_pos1, pat_pos2+1):
                if p in paternal_collapse[chrom].keys() and paternal_collapse[chrom][p]:
                    pat_tag += 1
            if pat_tag > 0:
                pat_sign = "minor"

            print(i + "\t" + mat_sign + "\t" + pat_sign)


if __name__ == '__main__':
    main()
