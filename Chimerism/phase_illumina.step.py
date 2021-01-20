#usr/bin/env python3
import sys

if len(sys.argv) < 4:
    sys.exit("python3 {} *.snp.reads win step".format(sys.argv[0]))

def snps2allele(snps):
    allele = ""
    pre = 0
    for a,b in snps:
        if pre == 0:
            allele = b
            pre = a
        else:
            if a == pre:
                # this is bad (a barcode support two genotypes in same position, e.g, 100,A; 100,T)
                return False
            else:
                allele += b
                pre = a

    return allele



def main():
    molecular = {}
    snps = set()
    records = {}
    # store all records in dict
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            line = line.strip()
            tmp = line.split("\t")
            pos = int(tmp[1])
            snps.add(pos)
            if pos not in records.keys():
                records[pos] = [line,]
            else:
                records[pos].append(line)

    # deal by windows and slide by step
    poses = sorted(list(snps))
    win = int(sys.argv[2])
    step = int(sys.argv[3])
    start = 0
    end = start + win
    while end <= len(poses):
        molecular = {}
        ref = []
        for p in poses[start:end]:
            for l in records[p]:
                tmp = l.split("\t")
                if (p, tmp[2]) not in ref:
                    ref.append((p, tmp[2]))
                if tmp[4] == 'ref':
                    base = tmp[2]
                else:
                    base = tmp[3]
                read_id = tmp[5]
                #barcode = read_id.split("_")[1]
                barcode = read_id
                if barcode not in molecular.keys():
                    molecular[barcode] = [(p, base),]
                else:
                    if (p, base) not in molecular[barcode]:
                        molecular[barcode].append((p, base))

        # different combinations of allele
        allele_types = {}
        for b in molecular.keys():
            if len(molecular[b]) == win:
                allele = snps2allele(molecular[b]) # link to a sequence
                if allele == False:
                    break
                # counting the sequence
                if allele not in allele_types.keys():
                    allele_types[allele] = 1
                else:
                    allele_types[allele] += 1

        allele_count = len(allele_types.keys())
        refseq = snps2allele(ref)
        if allele_count > 0:
            #print(">{}-{}\t{}\t{}".format(poses[start], poses[end-1], refseq, allele_count))
            pos_list = ",".join(map(lambda x: str(x), poses[start:end]))
            print(">{}\t{}\t{}".format(pos_list, refseq, allele_count))
            # sort alleles by frequency
            sorted_alleles = sorted(allele_types, key=lambda k: (allele_types[k], k), reverse=True)
            for k in sorted_alleles:
                print(k + "\t" + str(allele_types[k]))
        start += step
        end = start + win


if __name__ == '__main__':
    main()


