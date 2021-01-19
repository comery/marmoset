#!/usr/bin/env python3

import sys
if len(sys.argv) < 2:
    print("Usage: python3 {} {} ".format(sys.argv[0], "*.5evidence.xls"))
    exit()


def get_pattern(hap1, hap2):
    code = {'-': 0, 'A': 2, 'T': 3, 'G': 5, 'C': 7}
    return code[hap1] * code[hap2]



def call_haplotype(base_dis, ty=None):
    bases = {'A': int(base_dis[0]),
             'T': int(base_dis[1]),
             'G': int(base_dis[2]),
             'C': int(base_dis[3])}
    pi = float(base_dis[4])
    sorted_bases = sorted(bases, key=lambda k: bases[k], reverse=True)
    hap1 = sorted_bases[0]
    hap2 = sorted_bases[1]
    hap3 = sorted_bases[2]
    gt = {}
    if pi == 0:
        ##this is a homozygous site, bases[hap2] == 0
        return (hap1, bases[hap1], hap1, bases[hap1])
    elif ty == 'tenx':
        if pi <= 0.15: # this is a homozygous pattern
            if bases[hap2] >= 5 and bases[hap1] / bases[hap2] < 5:
                return False
            elif bases[hap2] >= 20:
                return (hap1, bases[hap1], hap2, bases[hap2])
            else:
                return (hap1, bases[hap1], hap1, bases[hap1]) #homozygous
        elif pi >= 0.4 and bases[hap3] >= 5:
            return False # what a mass!
        else:
            ## this is a heterozgyous pattern
            if bases[hap2] <= 3:
                return False
            else:
                return (hap1, bases[hap1], hap2, bases[hap2])
    else:
        return (hap1, bases[hap1], hap2, bases[hap2])


def which_wrong(tenx1, tenx2, mat, pat):
    # mat == pat
    if mat == pat:
        if mat == tenx1:
            return "eithor", tenx2
        elif mat == tenx2:
            return "eithor", tenx1
        else:
            return "both"
    else:
        if mat == tenx1 and pat != tenx2:
            return "pat", tenx2
        elif mat == tenx2 and pat != tenx1:
            return "pat", tenx1
        elif pat == tenx1 and mat != tenx2:
            return "mat", tenx2
        elif pat == tenx2 and mat != tenx1:
            return "mat", tenx1
        else:
            return 'both'


def backtracing(assembly_assumed, assembly, corrected, raw):
    cor_hap1, base_cor_hap1, cor_hap2, base_cor_hap2  = corrected
    raw_hap1, base_raw_hap1, raw_hap2, base_raw_hap2 = raw
    if assembly_assumed == assembly:
        """
        if assembly == cor_hap1:
            if assembly == raw_hap1:
                return "nothing"
            else:
                return "correcting_nice"
        else:
            if assembly == raw_hap1:
                return "polishing_nice"
            else:
                if cor_hap1 == raw_hap1:
                    return "polishing_nice"
                else:
                    return "other"
        """
        return "nothing"
    else:
        if assembly == cor_hap1:
            if cor_hap2 == assembly_assumed and base_cor_hap2 >= 5:
                return "assigning_error"
            else:
                # chain: assembly -> cor_hap1 -> raw_hap1
                return "sequencing_error" # error from original
        else:
            # chain: assembly ->x cor_hap1
            if cor_hap1 == assembly_assumed:
                if cor_hap2 == assembly:
                    return "assigning_error"
                else:
                    return "polishing_error"
            elif cor_hap2 == assembly_assumed and base_cor_hap2 >= 5:
                return "assigning_error"
            else:
                if cor_hap1 == raw_hap1 or cor_hap1 == raw_hap2:
                    return "sequencing_error"
                else:
                    return "complex_error"



class Record():
    def __init__(self, record):
        arr = record.strip().split("\t")
        self.pos = arr[1]
        self.matref = arr[2]
        self.patref = arr[4]
        self.mummer = arr[5]
        self.gatk = arr[6]
        self.note1 = arr[8]
        self.note2 = arr[9]
        self.tenx = arr[13:18]
        self.matCor = arr[21:26]
        self.matRaw = arr[29:34]
        self.patCor = arr[37:42]
        self.patRaw = arr[45:50]
        self.matcpos = arr[19]
        self.matrpos = arr[27]
        self.patcpos = arr[35]
        self.patrpos = arr[43]
        if arr[25] != "-" and arr[33] != "-" and arr[41] != "-" and arr[49] != "-":
            self.matpi = float(arr[25]) + float(arr[33])
            self.patpi = float(arr[41]) + float(arr[49])

    def get_tenx_info(self):
        return call_haplotype(self.tenx, ty='tenx')

    def get_matCor_info(self):
        return call_haplotype(self.matCor)

    def get_matRaw_info(self):
        return call_haplotype(self.matRaw)

    def get_patCor_info(self):
        return call_haplotype(self.patCor)

    def get_patRaw_info(self):
        return call_haplotype(self.patRaw)



def main():

    tag = open(sys.argv[1] + ".tag", 'w')
    with open(sys.argv[1], 'r') as fh:


        for i in fh:
            i = i.strip()
            if i.startswith("#"):
                tag.write(i + "\tINFO\tSNP\n")
                continue
            SNP = Record(i)
            
            if SNP.note1 == "GR" or SNP.note2 == "LOW":
                tag.write(i + "\thard2id\t-\n")
                continue


            # case 3, this site is located in grey region (GR) which is not on the
            # alignment of maternal/paternal genome, and 10x mapping in this site has
            # a low depth (<10), so it is hard to identity
            
            """
            if SNP.note2 == "LOW":
                # just see maternal evidence chain and paternal evidence, if both good
                if  SNP.matrpos == "-" or SNP.matcpos == "-" or SNP.patcpos == "-" or SNP.patrpos == "-":
                    tag.write(i + "\t" + "no_pacbio_read\t-\n")
                    continue
                error_type_mat = backtracing(SNP.matref, SNP.matref, SNP.get_matCor_info(), SNP.get_matRaw_info())
                error_type_pat = backtracing(SNP.patref, SNP.patref, SNP.get_patCor_info(), SNP.get_patRaw_info())
                if 'error' in error_type_mat and 'error' not in error_type_pat:
                    error_type = error_type_mat + "_mat" + "-" + "nothing_pat"
                elif 'error' in error_type_pat and 'error' not in error_type_mat:
                    error_type = "nothing_mat" + "-" + error_type_pat + "pat"
                else:
                    error_type = error_type_mat + "_mat" + "-" + error_type_pat + "_pat"

                if SNP.matref == SNP.patref:
                    tag.write(i + "\t" + error_type + "\t" + "False\n")
                else:
                    tag.write(i + "\t" + error_type + "\t" + "True\n")
                continue
            """
            # case 4, 
            if SNP.get_tenx_info() == False:
                tag.write(i + "\thard2id\t-\n")
                continue
            else:
                tenx_hap1, n1, tenx_hap2, n2 = SNP.get_tenx_info()


            # case 5, 
            """
            if SNP.note1 == "GR":
                # no paternal pacbio mapping result, so juct check maternal
                # backtracing(maternal) == 'nothing' and SNP.matref in [tenx1, tenx2]
                if  SNP.matrpos == "-" or SNP.matcpos == "-":
                    tag.write(i + "\t" + "no_pacbio_read\t-\n")
                    continue

                if SNP.matref in [tenx_hap1, tenx_hap2]:
                    error_type_mat = backtracing(SNP.matref, SNP.matref, SNP.get_matCor_info(), SNP.get_matRaw_info())
                    error_type = error_type_mat + "_mat"
                else:
                    error_type = "sequencing_error_mat"
                if tenx_hap1 == tenx_hap2:
                    tag.write(i + "\t" + error_type + "\t" + "False\n")
                else:
                    tag.write(i + "\t" + error_type + "\t" + "True\n")
                continue
            """

            
            if SNP.patref == SNP.gatk and SNP.mummer == "-" and tenx_hap1 != tenx_hap2:
                tag.write(i + "\tmummer_false_negative\tTrue\n")
                continue
            

            assembled_haplotypes = {SNP.matref, SNP.patref}
            tenx_haplotypes = {tenx_hap1, tenx_hap2}
            if assembled_haplotypes == tenx_haplotypes:
                tag.write(i + "\tNo_doubt\tTrue\n")
                continue

            """
            this an evidence chain to check which step error occur.
                                  10x pattern[h1,h2]
                                         ||
                    maternal      assembly pattern   paternal
                   /                                        \
                maternal corrected[h1,h2]  paternal corrected[h1,h2]
                /                                             \
            maternal raw[h1,h2]                     paternal raw[h1,h2]
            """
            # the assemblied haplotypes are not consistent with 10x haplotypes

            if which_wrong(tenx_hap1, tenx_hap2, SNP.matref, SNP.patref) == 'both':
                error_type = "terrible_error_both"
            else:
                # case 6, 
                if  SNP.matrpos == "-" or SNP.matcpos == "-" or SNP.patcpos == "-" or SNP.patrpos == "-":
                    tag.write(i + "\t" + "hard2id\t-\n")
                    continue

                who, shoulde_be = which_wrong(tenx_hap1, tenx_hap2, SNP.matref, SNP.patref)

                if who == 'mat':
                    # error occur in mat
                    error_type = backtracing(shoulde_be, SNP.matref, SNP.get_matCor_info(), SNP.get_matRaw_info())
                    error_type += "_mat"
                elif who == 'pat':
                    error_type = backtracing(shoulde_be, SNP.patref, SNP.get_patCor_info(), SNP.get_patRaw_info())
                    error_type += "_pat"
                else:
                    # eithor mat wrong or pat wrong, this is very complicated, you need check both
                    # of them, then determin which one is more likely wrong.

                    error_type_mat1 = backtracing(shoulde_be, SNP.matref, SNP.get_matCor_info(), SNP.get_matRaw_info())
                    error_type_pat1 = backtracing(SNP.patref, SNP.patref, SNP.get_patCor_info(), SNP.get_patRaw_info())

                    error_type_mat2 = backtracing(SNP.matref, SNP.matref, SNP.get_matCor_info(), SNP.get_matRaw_info())
                    error_type_pat2 = backtracing(shoulde_be, SNP.patref, SNP.get_patCor_info(), SNP.get_patRaw_info())


                    if 'error' in error_type_mat1 and 'error' not in error_type_pat1:
                        error_type_mat = backtracing(shoulde_be, SNP.matref,
                                                 SNP.get_matCor_info(),
                                                 SNP.get_matRaw_info())
                        error_type = error_type_mat + "_mat" + "-" + "nothing_pat"
                    elif 'error' not in error_type_mat2 and 'error' in error_type_pat2:
                        error_type_pat = backtracing(shoulde_be, SNP.patref,
                                                 SNP.get_patCor_info(),
                                                 SNP.get_patRaw_info())
                        error_type =  "nothing_mat" + "-" + error_type_pat + "_pat"
                    else:
                        error_type = "hard2id"
                #output
                if (tenx_hap1 == tenx_hap2):
                    tag.write(i + "\t" + error_type + "\t" + "False\n")
                else:
                    tag.write(i + "\t" + error_type + "\t" + "True\n")



if __name__ == '__main__':
    main()
