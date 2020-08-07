#!/usr/bin/env python3
import sys

if len(sys.argv) < 8:
    print("Usage: python3 {} {} {} {} {} {} {} {}".format(sys.argv[0],
                                                    "<compared.tab>",
                                                    "<rqpos.coords>",
                                                    "<10.base.xls>",
                                                    "<mat2mat.pb.cor.base.xls>",
                                                    "<mat2mat.pb.raw.base.xls>",
                                                    "<pat2pat.pb.cor.base.xls>",
                                                    "<pat2pat.pb.raw.base.xls>"))
    exit()

def get_base(file):
    data = {}
    with open(file ,'r') as fh:
        for i in fh.readlines():
            if i.startswith("Scaffold"):
                continue
            tmp = i.strip().split("\t")
            pos = tmp[1]
            info = tmp[0:3]
            info.extend(tmp[4:9])
            data[pos] = "\t".join(info)
    return data

def get_coords(file):
    coords = {}
    coords_all = {}
    with open(file ,'r') as fh:
        for i in fh.readlines():
            tmp = i.strip().split("\t")
            mat_pos = tmp[0]
            pat_pos = tmp[2]
            coords_all[mat_pos] = [tmp[2], tmp[3]]
            """
            some sites are not avaliably, because some site is on Grey Region(of not coordinative)
            and some sites are indels in Mummer alignment.
            """
            if pat_pos == '.' or pat_pos == "-":
                continue
            else:
                coords[mat_pos] = pat_pos
    return coords, coords_all

def get_info(pos, data, coords=None):
    non = "\t".join(['-'] * 8)
    if coords:
        # this is for paternal position, they based on maternal position
        if pos in coords.keys() and coords[pos] in data.keys():
            return data[coords[pos]]
        else:
            return non
    else:
        if pos in data.keys():
            return data[pos]
        else:
            return non

def main():
    coords, coords_all = get_coords(sys.argv[2])
    tenx = get_base(sys.argv[3])
    mat_cor_base = get_base(sys.argv[4])
    mat_raw_base = get_base(sys.argv[5])
    pat_cor_base = get_base(sys.argv[6])
    pat_raw_base = get_base(sys.argv[7])

    header = ["Chr", "Pos", "Depth", "A", "T", "G", "C", "Pi"]
    header = "\t".join(header * 5)
    tab_header = "#CHROME\tR_POS\tR_REF\tQ_POS\tQ_REF\tVAR_MUMMER\tVAR_GATK\tVAR_BCF\tNOTE_RG\tNOTE_DEP"
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#"):
                print(tab_header + "\t" + header)
                continue
            tab = line.split("\t")
            loca = tab[1]
            if loca in coords_all.keys():
                tab.insert(3, coords_all[loca][0])
                tab.insert(4, coords_all[loca][1])
            else:
                tab.insert(3, '-')
                tab.insert(4, '-')

            non = "\t".join(['-'] * 8)

            output = ["\t".join(tab),]
            output.append(get_info(loca, tenx))
            output.append(get_info(loca, mat_raw_base))
            output.append(get_info(loca, mat_cor_base))
            output.append(get_info(loca, pat_raw_base, coords))
            output.append(get_info(loca, pat_cor_base, coords))

            print("\t".join(output))

if __name__ == '__main__':
    main()

