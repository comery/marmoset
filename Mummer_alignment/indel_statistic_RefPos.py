import sys
import re
if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "out.indel", "prefix"))
    exit()

chrom = sys.argv[2]
out = open(chrom + ".indel.report.txt", 'w')

var_insert = {}
var_delete = {}
pos = []
kind = {}
pos_match = {}
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        i = i.strip()
        tmp = i.split("\t")
        if tmp[1] == ".":
            kind[tmp[0]] = "i"
            # inertion for Ref
            if tmp[0] not in var_insert.keys():
                pos.append(tmp[0])
                var_insert[tmp[0]] = tmp[2]
                pos_match[tmp[0]] = tmp[3]
            else:
                var_insert[tmp[0]] += tmp[2]

        elif tmp[2] == ".":
            kind[tmp[0]] = "d"
            # deletion for Ref
            if tmp[3] not in var_delete.keys():
                pos.append(tmp[0])
                var_delete[tmp[3]] = tmp[1]
                pos_match[tmp[0]] = tmp[3]
            else:
                var_delete[tmp[3]] += tmp[1]
        else:
            print("this is not a indel")

with open(sys.argv[2] + ".indel.report.txt", 'w') as out:
    print("#Chrom\tR\ttype\tsubType\tlength\tR_base\tQ_base", file=out)
    for p in pos:
        k = kind[p]
        if k == "i":
            r_pos = p
            q_pos = pos_match[p]
            indel_len = len(var_insert[r_pos])
            if indel_len == 1:
                t = "." + var_insert[r_pos]
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, r_pos, "ins", t,
                                                  1, ".", var_insert[r_pos]), file=out)
            else:
                q_pos1 = int(q_pos) + indel_len - 1
                t = ".S"
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, r_pos, "del", t,
                                                              indel_len, ".", var_insert[r_pos]), file=out)
        else:
            r_pos = p
            q_pos = pos_match[p]
            indel_len = len(var_delete[q_pos])
            if indel_len == 1:
                t = var_delete[q_pos] + "."
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, r_pos, "ins", t,
                                                  1, var_delete[q_pos], "."), file=out)
            else:
                r_pos1 = int(r_pos) + indel_len - 1
                t = "S."
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, r_pos, "del", t, indel_len,
                                                              var_delete[q_pos], "."), file=out)


