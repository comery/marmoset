#!/usr/bin/env python3
import os
import re
import sys
import subprocess
if len(sys.argv) < 5:
    note = """
    This script is to filter out these indels which don't have
    shotgun/pacbio reads supporting.
    """
    print(note)
    sys.exit("Usage: python3 {} {} {} {} {}".format(sys.argv[0],
                                              "indel.report",
                                              "*.base.xls",
                                              "output",
                                              "least_depth"))

least_depth = int(sys.argv[4])

def read_indel(f):
    ins = {}
    dels = {}
    indel_type = {}
    records = {}
    with open(f, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            else:
                tmp = i.strip().split("\t")
                if "-" in tmp[1]:
                    mat_pos = tmp[1].split("-")[0]
                    dels[int(mat_pos)] = tmp[6]
                elif "_" in tmp[1]:
                    mat_pos = tmp[1].split("_")[0]
                    ins[int(mat_pos)] = tmp[7]
                else:
                    mat_pos = tmp[1]
                    dels[int(mat_pos)] = tmp[6]

                indel_type[int(mat_pos)] = tmp[3]
                records[int(mat_pos)] = i

    return ins, dels, records


def is_support(indel, tenx):
    sign = 0
    for t in tenx.split(";"):
        tenx_pattern= re.compile(r'(\S+)\(')
        tenx_indel = tenx_pattern.findall(t)
        number_pattern = re.compile(r'\S+(\d+)')
        number = number_pattern.findall(t)
        if indel == tenx_indel[0] and int(number[0]) >= least_depth:
            sign = 1
    if sign:
        return True
    else:
        return False

def main():
    ins, dels, records = read_indel(sys.argv[1])
    unmapped_indel = 0
    with open(sys.argv[2], 'r') as fh, open(sys.argv[3], 'w') as fw:
        for r in fh:
            if r.startswith("Scaffold"):
                continue
            if len(ins.keys()) == 0 and len(dels.keys()) == 0:
                break
            tmp = r.strip().split("\t")
            mat_pos = int(tmp[1])
            real_dep = int(tmp[11])
            insert = tmp[9]
            delete = tmp[10]
            pos_for_del = mat_pos + 1
            if mat_pos in ins.keys() or mat_pos in dels.keys():
                if real_dep < 3:
                    unmapped_indel += 1
                    continue
            if insert and mat_pos in ins.keys():
                if is_support(ins[mat_pos], insert) == True:
                    fw.write(records[mat_pos])
                del ins[mat_pos]
            elif delete and pos_for_del in dels.keys():
                if is_support(dels[pos_for_del], delete) == True:
                    fw.write(records[pos_for_del])
                del dels[pos_for_del]

    sys.stdout.write("unmapped_indel:\t{}".format(unmapped_indel))
if __name__ == '__main__':
    main()
