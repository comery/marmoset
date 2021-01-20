#usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit("python3 {} chr1.genotype.all.reads.4_1.alleles".format(sys.argv[0]))



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
        minior_region = 0
        sequencing_err_region = 0
        total_region = 0
        minior_snps = 0
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split("\t")
            poses = tmp[0].replace(">", "").split("-")
            poses = list(map(lambda x: int(x), poses))
            total_region += poses[-1] - poses[0] + 1
            if tmp[-1] == "minor":
                minior_region += poses[-1] - poses[0] + 1
                minior_snps += len(poses)
            elif tmp[-1] == "sequencing_err":
                sequencing_err_region += poses[-1] - poses[0] + 1
        sequencing_err_rate = sequencing_err_region / total_region
        minior_rate = minior_region / total_region
        print("#total_region\tsequencing_err_region\tsequencing_err_rate\tminior_region\tminior_rate\tminior_snps")
        print("{}\t{}\t{:.6f}\t{}\t{:.6f}\t{}".format(total_region, sequencing_err_region,
                                              sequencing_err_rate,
                                              minior_region,
                                              minior_rate,
                                              minior_snps))



if __name__ == '__main__':
    main()
