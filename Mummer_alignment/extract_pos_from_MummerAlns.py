#!/usr/bin/env python3
import sys

usage = "this script is to extract ref/query coordinative position, you'd better\n" +\
        "use one to one coords alignment which generated by *.1coords file.\n"
if len(sys.argv) < 4:
    print(usage)
    print("Usage: python3 {} {} {} {}".format(sys.argv[0], "<pos>", "<mummer.alns>", "<outpre>"))
    exit()

def read_pos(file):
    pos = {}
    with open(file, 'r') as fh:
        for p in fh.readlines():
            pos[int(p.strip())] = 1
    return pos


def get_info(r_pos, r_seq, q_pos, q_seq):
    rbase = {}
    qbase = {}
    qdels = {}
    rqcoords = {}
    if len(r_seq) != len(q_seq):
        print("bad line with different length")
        exit()
    r_skip = 0
    q_skip = 0
    for i in range(len(r_seq)):
        rb = r_seq[i]
        qb = q_seq[i]

        if rb == ".":
            # rb is deletion, qb is 'a, t, c, g n';
            r_skip += 1

        elif qb == ".":
            # rb is 'a, t, c, g n'; qb is deletion
            q_skip += 1
            r_real_pos = r_pos + i  - r_skip
            qdels[r_real_pos] = rb.upper()

        else:
            # rb, qb is 'a, t, c, g n'
            r_real_pos = r_pos + i  - r_skip
            rbase[r_real_pos] = rb.upper()
            q_real_pos = q_pos + i - q_skip
            qbase[q_real_pos] = qb.upper()
            rqcoords[r_real_pos] = q_real_pos

    return (rbase, qbase, qdels, rqcoords)

def main():
    allrecords = {}
    snps = read_pos(sys.argv[1])
    pass_symbol = [' ', '-', '=', '/']
    with open(sys.argv[2], 'r') as fh:
        for line in fh:
            if line[0] in pass_symbol:
                continue
            elif '^' in line:
                continue
            else:
                line = line.strip()
            tmp = line.split()
            if len(tmp) == 2:
                line1 = fh.readline().strip()
                tmp1 = line1.split()
                r_pos = int(tmp[0])
                r_seq = tmp[1]
                q_pos = int(tmp1[0])
                q_seq = tmp1[1]
                (rbase, qbase, qdels, rqcoords) = get_info(r_pos, r_seq, q_pos, q_seq)
                allsites = list(rbase.keys()) + list(qdels.keys())
                for b in allsites:
                    if b in snps.keys():
                        if b in qdels.keys():
                            allrecords[b] = "{}\t{}\t{}\t{}".format(b, qdels[b], ".", ".")
                            del snps[b]
                        else:
                            allrecords[b] = "{}\t{}\t{}\t{}".format(b, rbase[b], rqcoords[b],
                                                                qbase[rqcoords[b]])
                            del snps[b]
        for e in snps.keys():
            allrecords[e] = str(e) + "\t-\t-\t-\t"

    with open(sys.argv[3] + ".rqpos.coords", 'w') as fh:
        for k in sorted(allrecords.keys()):
            print(allrecords[k], file=fh)


if __name__ == '__main__':
    main()
