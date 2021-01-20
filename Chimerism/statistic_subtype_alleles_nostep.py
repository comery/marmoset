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
        phased_region_len = 0
        one_alleles_region = 0
        two_alleles_region = 0
        three_alleles_region = 0
        four_alleles_region = 0
        other_alleles_region = 0
        for i in fh:
            if not i.startswith(">"):
                continue
            tmp = i.strip().split("\t")
            poses = tmp[0].replace(">", "").split(",")
            poses = list(map(lambda x: int(x), poses))
            refseq = tmp[1]
            type_number = int(tmp[2])
            phased_region_len += poses[-1] - poses[0] + 1
            if type_number == 1:
                one_alleles_region += poses[-1] - poses[0] + 1
            if type_number == 2:
                two_alleles_region += poses[-1] - poses[0] + 1
            if type_number == 3:
                three_alleles_region += poses[-1] - poses[0] + 1
            elif type_number == 4:
                four_alleles_region += poses[-1] - poses[0] + 1
            elif type_number > 4:
                other_alleles_region += poses[-1] - poses[0] + 1

        title = ['phased_region_len', 'one_alleles_region', 'one_alleles_rate', 'two_alleles_region', 'two_alleles_rate', 'three_alleles_region', 'three_type_rate', 'four_alleles_region', 'four_alleles_rate', 'other_alleles_region']
        print("#" + "\t".join(title))
        print("{}\t{}\t{:.4f}\t{}\t{:.4f}\t{}\t{:.4f}\t{}\t{:.4f}\t{}\t{:.4f}".format(phased_region_len,
                                                                  one_alleles_region,
                                                                  one_alleles_region / phased_region_len,
                                                                  two_alleles_region,
                                                                  two_alleles_region / phased_region_len,
                                                                  three_alleles_region,
                                                                  three_alleles_region / phased_region_len,
                                                                  four_alleles_region,
                                                                  four_alleles_region / phased_region_len,
                                                                  other_alleles_region,
                                                                  other_alleles_region / phased_region_len))


if __name__ == '__main__':
    main()
