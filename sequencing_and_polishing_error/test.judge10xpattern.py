import sys
def call_haplotype(base_dis, ty=None):
    bases = {'A': int(base_dis[0]),
             'T': int(base_dis[1]),
             'G': int(base_dis[2]),
             'C': int(base_dis[3])}
    pi = float(base_dis[4])
    sorted_bases = sorted(bases, key=lambda k: bases[k], reverse=True)
    hap1 = sorted_bases[0]
    hap2 = sorted_bases[1]
    if pi == 0:
        hap2 = hap1 ##this is a homozygous site, bases[hap2] == 0
        return hap1, hap2
    elif ty:
        if pi < 0.15: # this is a homozygous pattern
            if bases[hap2] >= 5 and bases[hap1] / bases[hap2] < 5:
                return False
            else:
                return hap1, hap1 #homozygous
        else:
            ## this is a heterozgyous pattern
            if bases[hap2] <= 3:
                return False
            else:
                return hap1, hap2
    else:
        return hap1, hap2

def main():
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            i = i.strip()
            if i.startswith("A"):
                print(i)
                continue
            elif '-' in i:
                continue
            tmp = i.split()
            if call_haplotype(tmp, 'tenx') == False:
                print(i + "\thard2id")
            else:
                hap1, hap2 = call_haplotype(tmp, 'tenx')
                print(i + "\t" + hap1 + "\t" + hap2)


if __name__ == '__main__':
    main()
