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

def main():
    with open(sys.argv[1], 'r') as fh:
        print("#block\tsnps\tallele_type\tref_allele\tref_count\tref_rate\tref_corrected_rate\tnote")
        for i in read_fasta(fh):
            name, stat = i
            name = name.replace(">", "")
            result = ""
            tmp = name.split("\t")
            poses = tmp[0].split(",")
            poses_info = "(" + ",".join(poses) + ")"
            ref_allele = tmp[1]
            allele_type = int(tmp[2])
            potential_alleles = {}
            for s in stat:
                tmp1 = s.split("\t")
                potential_alleles[tmp1[0]] = int(tmp1[1].strip())
            sorted_alleles = sorted(potential_alleles, key=lambda k: (potential_alleles[k], k), reverse=True)
            """
            if ref_allele in sorted_alleles[0:2] and potential_alleles[1] > potential_alleles[2]:
                result = "major"
            elif ref_allele in sorted_alleles[0:3] and potential_alleles[1] == potential_alleles[2]:
                result = "major"
            elif ref_allele in sorted_alleles[0:3] and potential_alleles[0] == potential_alleles[1] and potential_alleles[1] == potential_alleles[2] and potential_alleles[2] > potential_alleles[3]:
                result = "major"
            elif ref_allele in sorted_alleles[0:4] and potential_alleles[0] == potential_alleles[3] and potential_alleles[1] > potential_alleles[2]:
                result = "unclear"
            else:
                result = "chimeric"
            """
            ref_rate = get_rate(ref_allele, potential_alleles)
            ref_corrected_rate = 0
            if ref_allele in potential_alleles.keys():
                ref_count = potential_alleles[ref_allele]
            else:
                ref_count = 0
            allele_number = len(sorted_alleles)
            if ref_allele not in sorted_alleles:
                result = "sequencing_err"
            elif allele_number == 1:
                result = "major"
                ref_corrected_rate = ref_rate
            elif allele_number == 2:
                result = "major"
                ref_corrected_rate = ref_rate
            elif allele_number == 3:
                if potential_alleles[sorted_alleles[-1]] > 1:
                    if ref_rate >= 0.333 or ref_allele in sorted_alleles[0:2]:
                        result = "major"
                        ref_corrected_rate = ref_rate
                    else:
                        result = "minor"
                        ref_corrected_rate = ref_rate
                elif potential_alleles[sorted_alleles[-1]] == 1:
                    if ref_allele == sorted_alleles[-1]:
                        result = "sequencing_err"
                    else:
                        result = "major"
                        ref_corrected_rate = ref_rate
            elif allele_number > 3:
                #filer allele of number = 1
                new_alleles = list(filter(lambda x:potential_alleles[x] > 1, sorted_alleles))
                if len(new_alleles) > 0:
                    cutoff = 1 / len(new_alleles)
                    total_filtered = 0
                    for k in new_alleles:
                        total_filtered += potential_alleles[k]
                    ref_corrected_rate = potential_alleles[ref_allele] / total_filtered
                    if ref_allele not in new_alleles:
                        result = "sequencing_err"
                    elif ref_corrected_rate >= cutoff:
                        result = "major"
                    else:
                        result = "minor"

                else:
                    result = "unclear"
            print("{}-{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{}".format(poses[0], poses[-1], poses_info, allele_type, ref_allele, ref_count, ref_rate, ref_corrected_rate, result))



if __name__ == '__main__':
    main()
