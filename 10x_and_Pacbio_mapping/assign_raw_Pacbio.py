#!/usr/bin/env python3
import os
import sys
import gzip

if len(sys.argv) < 5:
    print("Usage: python3 {} {} {} {} {}".format(sys.argv[0] , "raw.fasta", "mat.id", "pat.id", "outdir"))
    exit()

def get_id(file):
    data = {}
    with open(file, 'r') as fh:
        for i in fh:
            record = i.strip().split()[0]
            data[record] = 1
    return data

def addtwodimdict(thedict, key_a, key_b):
    if key_a in thedict:
        thedict[key_a].update({key_b: thedict[key_a][key_b] + 1})
    else:
        thedict.update({key_a:{key_b: 1}})

def main():
    statis = {}
    items = ['mat', 'pat', 'both', 'none']
    mat_ids = get_id(sys.argv[2])
    pat_ids = get_id(sys.argv[3])
    outdir = sys.argv[4]
    from Bio import SeqIO

    rawfasta = sys.argv[1].strip()
    if rawfasta.endswith(".gz"):
        out = gzip.open(rawfasta, 'rt')
    else:
        out = open(rawfasta, 'r')
    filename = rawfasta.split("/")[-1]
    basename = filename.replace(".subreads.bam.fasta.gz", "")
    #print(filename + "\t" + basename)
    statis[filename] = {}
    for t in items:
        statis[filename].update({t: 0})
    mat_assigned = outdir + "/mat/" + basename + ".fa.gz"
    pat_assigned = outdir + "/pat/" + basename + ".fa.gz"
    with gzip.open(mat_assigned, 'wt') as mh, gzip.open(pat_assigned, 'wt') as ph:
        for record in SeqIO.parse(out, "fasta"):
            name = str(record.id)
            seq = str(record.seq)
            if name in mat_ids.keys() and name not in pat_ids.keys():
                mh.write(">" + name + "\n" + seq + "\n")
                addtwodimdict(statis, filename, 'mat')
                #statis[filename]['mat'] += 1
            elif name in pat_ids.keys() and name not in mat_ids.keys():
                ph.write(">" + name + "\n" + seq + "\n")
                addtwodimdict(statis, filename, 'pat')
                #statis[filename]['pat'] += 1
            elif name in mat_ids.keys() and name in pat_ids.keys():
                #statis[filename]['both'] += 1
                addtwodimdict(statis, filename, 'both')
            else:
                #statis[filename]['none'] += 1
                addtwodimdict(statis, filename, 'none')


    out = [filename,]
    for i in items:
        out.append(str(statis[filename][i]))
    with open(outdir + "/" + basename + ".log", 'w') as fw:
        print("\t".join(out), file=fw)

if __name__ == '__main__':
    if os.path.exists(sys.argv[4]) == False:
        print("outdir does not exists")
        exit()
    elif os.path.exists(sys.argv[4] + "/mat") == False:
        print("outdir {}/mat does not exists".format(sys.argv[4]))
        exit()
    elif os.path.exists(sys.argv[4] + "/pat") == False:
        print("outdir {}/pat does not exists".format(sys.argv[4]))
        exit()
    main()

