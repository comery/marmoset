#!/usr/bin/env python3

import sys
if len(sys.argv) < 2:
    print("Usage: python3 {} {} ".format(sys.argv[0], "*.base.xls.check"))
    exit()

def get_info(arr):
    """
    #CHROME    R_POS   R_REF   Q_POS   Q_REF   VAR_MUMMER  VAR_GATK    VAR_BCF NOTE_RG NOTE_DEP
    Chr Pos Depth   A   T   G   C   Pi  Chr Pos Depth   A   T   G   C   Pi  Chr Pos Depth   A   T
    G   C   Pi
    """
    matref = arr[2]
    patref = arr[4]
    tenx_bases = arr[13:18]
    mat_bases = arr[21:26]
    pat_bases = arr[-5:]

    return (matref, patref, tenx_bases, mat_bases, pat_bases)

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
        ##this is a homozygous site, bases[hap2] == 0
        return hap1, hap1
    elif ty == 'tenx':
        #print("{}:{}\t{}:{}".format(hap1, bases[hap1], hap2, bases[hap2]))
        if pi <= 0.15: # this is a homozygous pattern
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
    irint(base_dis)

def main():
    snp = open("SNP.checked.true.xls", 'w')
    tag = open(sys.argv[1] + ".tag", 'w')
    with open(sys.argv[1], 'r') as fh:
        code = {'A': 2, 'T': 3, 'G': 5, 'C': 7}

        # init all type
        all_potential = 0
        grey_region = 0
        low_depth = 0
        no_doubt_good = 0
        false_indel = 0
        mummer_false_negative = 0
        no_pacbio_read = 0
        hard_to_identify = 0
        gatk_false_negative = 0
        sequencing_err = 0
        polishing_err = 0
        complex_err = 0
        other = 0


        for i in fh:
            i = i.strip()
            if i.startswith("#"):
                snp.write(i + "\n")
                tag.write(i + "\tINFO\n")
                continue
            arr = i.split("\t")
            # basic static
            (matref, patref, tenx_bases, mat_bases, pat_bases) = get_info(arr)


            mummer_var = arr[5]
            gatk_var = arr[6]
            sambcf_var = arr[7]
            note_gr = arr[8]
            note_dep = arr[9]

            all_potential += 1
            if note_gr == "GR":
                grey_region += 1
                tag.write(i + "\tGR\n")
                continue
            if note_dep == "LOW":
                low_depth += 1
                tag.write(i + "\tLOW\n")
                continue
            if patref == ".":
                false_indel += 1
                tag.write(i + "\tFalse_Indel_True_Snp\n")
                snp.write(i + "\n")
                continue
            if patref == mummer_var and mummer_var == gatk_var and gatk_var == sambcf_var:
                no_doubt_good += 1
                tag.write(i + "\tNo_doubt\n")
                snp.write(i + "\n")
                continue
            # mat2mat corrected pb reads;
            # pat2pat corrected pb reads;
            if arr[20] == "-" or arr[-6] == "-":
                tag.write(i + "\t" + "no_pacbio_read\n")
                no_pacbio_read += 1
                continue
            #print("assembly: {} {}".format(matref, patref))

            if call_haplotype(tenx_bases, ty='tenx') == False:
                tag.write(i + "\t" + "hard_to_identify\n")
                hard_to_identify += 1
                continue
            else:
                tenx_hap1, tenx_hap2 = call_haplotype(tenx_bases, ty='tenx')

            if patref == gatk_var and mummer_var == "-" and tenx_hap1 != tenx_hap2:
                mummer_false_negative += 1
                tag.write(i + "\tmummer_false_negative\n")
                snp.write(i + "\n")
                continue

            if tenx_hap1 != tenx_hap2:
                # heterogyous pattern based on 10x data
                snp.write(i + "\n")

            assembly_pattern = code[matref] * code[patref]
            tenx_pattern = code[tenx_hap1] * code[tenx_hap2]
            mat_hap_cor, na = call_haplotype(mat_bases)
            pat_hap_cor, na = call_haplotype(pat_bases)
            if assembly_pattern == tenx_pattern:
                if (gatk_var == "-" and matref == mat_hap_cor
                    and patref == pat_hap_cor and tenx_hap1!= tenx_hap2):
                    tag.write(i + "\t" + "gatk_false_negative\n")
                    #snp.write(i + "\n") # included above
                    gatk_false_negative += 1

                else:
                    tag.write(i + "\t" + "complex_err\n")
                    complex_err += 1
            elif matref == tenx_hap1 and patref != tenx_hap2 and patref == pat_hap_cor:
                sequencing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "pat_sequencing_err\n")
                else:
                    tag.write(i + "\t" + "sequencing_err\n")
            elif matref == tenx_hap1 and patref != tenx_hap2 and patref != pat_hap_cor:
                polishing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "pat_polishing_err\n")
                else:
                    tag.write(i + "\t" + "polishing_err\n")
            elif matref == tenx_hap2 and patref != tenx_hap1 and patref == pat_hap_cor:
                sequencing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "pat_sequencing_err\n")
                else:
                    tag.write(i + "\t" + "sequencing_err\n")
            elif matref == tenx_hap2 and patref != tenx_hap1 and patref != pat_hap_cor:
                polishing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "pat_polishing_err\n")
                else:
                    tag.write(i + "\t" + "polishing_err\n")

            elif matref != tenx_hap1 and patref == tenx_hap2 and matref == mat_hap_cor:
                sequencing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "mat_sequencing_err\n")
                else:
                    tag.write(i + "\t" + "sequencing_err\n")
            elif matref != tenx_hap1 and patref == tenx_hap2 and matref != mat_hap_cor:
                polishing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "mat_polishing_err\n")
                else:
                    tag.write(i + "\t" + "polishing_err\n")
            elif matref != tenx_hap2 and patref == tenx_hap1 and matref == mat_hap_cor:
                sequencing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "mat_sequencing_err\n")
                else:
                    tag.write(i + "\t" + "sequencing_err\n")
            elif matref != tenx_hap2 and patref == tenx_hap1 and matref != mat_hap_cor:
                polishing_err += 1
                if matref != patref:
                    tag.write(i + "\t" + "mat_polishing_err\n")
                else:
                    tag.write(i + "\t" + "polishing_err\n")
            elif (matref != tenx_hap1 and patref != tenx_hap2) and (matref != tenx_hap2 and patref !=tenx_hap1):
                tag.write(i + "\t" + "complex_err\n")
                complex_err += 1
            else:
                tag.write(i + "\t" + "other\n")
                other += 1

    snp.close()
    tag.close()
    print("all_potential: {}".format(all_potential))
    print("grey_region: {}".format(grey_region))
    print("low_depth: {}".format(low_depth))
    print("no_doubt_good: {}".format(no_doubt_good))
    print("False_Indel_True_Snp: {}".format(false_indel))
    print("mummer_false_negative: {}".format(mummer_false_negative))
    print("no_pacbio_read: {}".format(no_pacbio_read))
    print("hard_to_identify: {}".format(hard_to_identify))
    print("gatk_false_negative: {}".format(gatk_false_negative))
    print("sequencing_err: {}".format(sequencing_err))
    print("polishing_err: {}".format(polishing_err))
    print("complex_err: {}".format(complex_err))
    print("other: {}".format(other))


if __name__ == '__main__':
    main()
