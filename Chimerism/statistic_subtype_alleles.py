#usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit("python3 {} chr1.genotype.all.reads.4_1.alleles".format(sys.argv[0]))


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, seq)

def get_rate(seq, alleles):
    total_numer = 0
    for k in alleles.keys():
        total_numer += alleles[k]
    if seq in alleles.keys():
        rate = alleles[seq] / total_numer
    else:
        rate = 0

    return rate

def overlaped(list1, list2):
    if list2[0] <= list1[0]:
        return True
    else:
        return False

def main():
    with open(sys.argv[1], 'r') as fh:
        pre = []
        for i in fh:
            tmp = i.strip().split("\t")
            poses = tmp[0].replace(">", "").split(",")
            refseq = tmp[1]
            type_number = int(tmp[2])
            phased_region_len = 0
            three_alleles_region = 0
            four_alleles_region = 0
            mix_three_four_region = 0
            other_alleles_region = 0
            tmp = []
            tmp_allele_types = []
            if not pre:
                #phased_region_len = poses[-1] - poses[0] + 1
                tmp = poses
                tmp_allele_types.append(type_number)
                pre = poses
            else:
                if overlaped(pre, poses):
                    tmp = pre.append(poses[-1]) # if overlaped, the second has one more snp (-1) because step=1;
                else:
                    #output tmp, classfy




            phased_region_len = solve_overlap(phased_region_len, poses)




if __name__ == '__main__':
    main()
