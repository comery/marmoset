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

def get_min_depth(stat):
    tmp_depth = []
    for s in stat:
        tmp = s.split("\t")
        tmp_depth.append(int(tmp[-1]))
    tmp_depth.sort()
    if tmp_depth[0] == 1 and tmp_depth[1] > 1:
        return 1
    else:
        return 1.5


def main():
    with open(sys.argv[1], 'r') as fh:
        (count_type3, len_type3, count_type3_err, len_type3_err) = (0, 0, 0, 0)
        (count_type4, len_type4, count_type4_err, len_type4_err) = (0, 0, 0, 0)

        for i in read_fasta(fh):
            name, stat = i
            name = name.replace(">", "")
            result = ""
            tmp = name.split("\t")
            poses = tmp[0].split(",")
            win_len = int(poses[-1]) - int(poses[0]) + 1
            ref_allele = tmp[1]
            allele_type = int(tmp[2])
            if allele_type == 3:
                tmp_depth = []
                len_type3 += win_len
                count_type3 += 1
                if get_min_depth(stat) == 1:
                    len_type3_err += win_len
                    count_type3_err += 1
                    print(">" + name + "\n" + "\n".join(stat))


            elif allele_type == 4:
                tmp_depth = []
                len_type4 += win_len
                count_type4 += 1
                if get_min_depth(stat) == 1:
                    len_type4_err += win_len
                    count_type4_err += 1
                    print(">" + name + "\n" + "\n".join(stat))
        print("3:\t{}\t{}\t{}\t{}".format(count_type3, len_type3, count_type3_err, len_type3_err))
        print("4:\t{}\t{}\t{}\t{}".format(count_type4, len_type4, count_type4_err, len_type4_err))
if __name__ == '__main__':
    main()
