#!/usr/bin/env python3
import re
import sys
import argparse

usage = "this script is to design primer for SNP validation, acorrding\n" +\
        "to SNP location, I will extract up/down stream region to design\n" +\
        "primers\n" +\
        "check primer follow behind standards and priorities:\n" +\
        "- containing target SNP site\n" +\
        "- no flanking indel or other SNP\n" +\
        "- ploy N in primers, and no GGG/CCC/A on edge of 3' end\n" +\
        "- avoid GC/AT rich in 3' end" +\
        "- bad Tm > max_tm or < min_tm or > tm_differene\n" +\
        "- hairpin structure\n" +\
        "- dimer structure\n"

parser = argparse.ArgumentParser(description=usage,
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-in", metavar='<File>', type=str, dest='infile', required=True,
                    help='input file, Mummer alns or aligned fasta')
parser.add_argument("-pos", metavar='<File>', type=str, dest='position', required=True,
                    help='target position you want to amplify')
parser.add_argument("-out", metavar='<STR>', type=str, dest='outpre', default='out',
                    help='out prefix for output files')
parser.add_argument("-min_size", metavar='<INT>', type=int, dest='min_size',
                    help='min amplicon length you expected, default=400', default=400)
parser.add_argument("-max_size", metavar='<INT>', type=int, dest='max_size',
                    help='max amplicon length you expected, default=600', default=600)
parser.add_argument("-min_plen", metavar='<INT>', type=int, dest='min_plen',
                    help='min primer length you expected, default=20', default=20)
parser.add_argument("-max_plen", metavar='<INT>', type=int, dest='max_plen',
                    help='max primer length you expected, default=25', default=25)
parser.add_argument("-min_gc", metavar='<Float>', type=float, dest='min_gc',
                    help='min GC content of primer, default=0.4', default=0.4)
parser.add_argument("-max_gc", metavar='<Float>', type=float, dest='max_gc',
                    help='max GC content of primer, default=0.6', default=0.6)
parser.add_argument("-min_tm", metavar='<INT>', type=int, dest='min_tm',
                    help='amplicon length you expected, default=55', default=55)
parser.add_argument("-max_tm", metavar='<INT>', type=int, dest='max_tm',
                    help='amplicon length you expected, default=75', default=75)
parser.add_argument("-tm_diff", metavar='<INT>', type=int, dest='tm_difference',
                    help='tm difference between forward and reverse primer, default=4', default=4)
parser.add_argument("-max_continous_dimer", metavar='<INT>', type=int, dest='max_continous_dimer',
                    help='max continous base for detecting dimer/hairpin, default=4', default=4)
parser.add_argument("-extra_snp_allow", metavar='<INT>', type=int, dest='extra_snp_allow',
                    help='extra snp allowed in applicon region expect target SNP, default=0', default=0)

parser.add_argument("-step", metavar='<INT>', type=int, dest='step',
                   help='step length when scaning whole alignment, default=1', default=1)


def even2odd(number):
    if number % 2 == 0:
        return number + 1
    else:
        return number

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()
if args.max_continous_dimer < 3:
    print("-max_continous_dimer={} is too strict, please change!".format(args.max_continous_dimer))
    exit()

if args.min_size >= args.max_size:
    print("-min_size={} and -max_size={} are invalid, please cheange".format(args.min_size, args.max_size))
    exit()
else:
    args.min_size = even2odd(args.min_size)
    args.max_size = even2odd(args.max_size)

if args.max_size < args.min_plen * 2:
    print("too short to design!")
    exit()



def get_phase(string):
    # string is like: --   END alignment [ +1 5334 - 61486 | +1 4 - 56173 ]
    # and get numer '5334' and '4'
    pattern_start = re.compile(r'(\d+) -')
    pattern_end = re.compile(r'- (\d+)')
    phase = {}
    end = {}
    m0 = pattern_start.findall(string)
    m1 = pattern_end.findall(string)
    if m0:
        phase['mat'] = int(m0[0])
        phase['pat'] = int(m0[1])
    if m1:
        end['mat'] = int(m1[0])
        end['pat'] = int(m1[1])

    return phase, end

def get_alignment(fh):
    def __alignment_end__(string):
        pattern = re.compile(r'END')
        if pattern.search(string):
            return True
        else:
            return False

    align = {'mat': '', 'pat': ''}
    pass_symbol = ['=', '/']

    for line in fh:
        if line[0] in pass_symbol or len(line.strip()) == 0:
            continue
        elif '^' in line or 'Alignments' in line or 'BEGIN' in line:
            continue
        else:
            line = line.strip()
            tmp = line.split()
            if __alignment_end__(line):
                phase, end = get_phase(line)
                if align:
                    yield align, phase, end
                    align = {'mat': '', 'pat': ''}
            elif len(tmp) == 2:
                # line is mat pos and seq
                # line1 is pat pos and seq
                line1 = fh.readline().strip()
                tmp1 = line1.split()
                r_seq = tmp[1]
                q_seq = tmp1[1]
                align['mat'] += r_seq
                align['pat'] += q_seq

    return align, phase, end  # last item


def slide_win(aln, phase, win_len, step=args.step):
    mat_pos = phase['mat']
    pat_pos = phase['pat']
    mat_seq = aln['mat']
    pat_seq = aln['pat']
    skip_mat = 0
    skip_pat = 0
    i = 0
    while i < len(mat_seq) - win_len:
        tmp_mat_seq = mat_seq[i:i+win_len]
        tmp_pat_seq = pat_seq[i:i+win_len]
        if i == 0:
            real_start_pos_mat = mat_pos
            real_start_pos_pat = pat_pos
        else:
            walked_mat = mat_seq[i-step:i]
            walked_pat = pat_seq[i-step:i]
            skip_m = walked_mat.count(".")
            skip_p = walked_pat.count(".")
            skip_mat += skip_m
            skip_pat += skip_p
            real_start_pos_mat = mat_pos - skip_mat + i
            real_start_pos_pat = pat_pos - skip_pat + i
        tmp_phase = {'mat': real_start_pos_mat, 'pat': real_start_pos_pat}
        tmp_aln = {'mat': tmp_mat_seq, 'pat': tmp_pat_seq}
        yield tmp_aln, tmp_phase
        i += step


def matched(s1, s2, length):
    code = {'A':'T',
            'T':'A',
            'C':'G',
            'G':'C',
            '-': 'N'}
    continuous_matched = 0
    all_continuous_matches = [0,]
    for i in range(length):
        if s1[i] == code[s2[i]]:
            if continuous_matched == 0:
                continuous_matched = 1
            else:
                continuous_matched += 1
        else:
            all_continuous_matches.append(continuous_matched)
            continuous_matched = 0
    longest_continuous_matched = max(all_continuous_matches)
    return longest_continuous_matched

def reverse_complementation(sequence):
    transtable = str.maketrans('ATGC-', 'TACG-')
    sequence = sequence.translate(transtable)
    return sequence[::-1]

def complementation(sequence):
    transtable = str.maketrans('ATGC-', 'TACG-')
    sequence = sequence.translate(transtable)
    return sequence

def ploy_x(seq):
    pa = re.compile(r'A{4,}')
    pt = re.compile(r'T{4,}')
    pc = re.compile(r'C{4,}')
    pg = re.compile(r'G{4,}')
    pat = "AT" * 4
    pac = "AC" * 4
    pag = "AG" * 4
    ptc = "TC" * 4
    ptg = "TG" * 4
    pcg = "CG" * 4
    pta = "TA" * 4
    pca = "CA" * 4
    pga = "GA" * 4
    pct = "CT" * 4
    pgt = "GT" * 4
    pgc = "GC" * 4
    p3endG = re.compile(r'GGG$')
    p3endC = re.compile(r'CCC$')
    #p3endA = re.compile(r'A$')
    if pa.search(seq) or pt.search(seq) or pc.search(seq) or pg.search(seq):
        return True
    elif (pat in seq or pac in seq or pag in seq or ptc in seq or ptg in seq or pcg in seq
        or pta in seq or pca in seq or pga in seq or pct in seq or pgt in seq or pgc in
          seq):
        return True
    elif p3endG.search(seq) or p3endC.search(seq):
        return True
    else:
        return False



def dimer(primer1, primer2, continuous_matched_allow=args.max_continous_dimer):
    # chech how many matched bases(complement) between two sequences

    # repair short sequence to have same length with longer one
    if len(primer1) >= len(primer2):
        repair = (len(primer1) - len(primer2)) * '-'
        longer_len = len(primer1)
        primer2 += repair
    else:
        repair = (len(primer2) - len(primer1)) * '-'
        longer_len = len(primer2)
        primer1 += repair

    """
    scanning
    ----->5'-3'
    atcgagaggact
             |||
             tgacgaggggat
                 3'<----5'
    ----->5'-3'
    atcgagaggact
        |     |
        tgacgaggggat
            3'<----5'
    """
    signal = 0
    for i in range(3,longer_len):
        s1 = primer1[-i:]
        s2 = primer2[0:i]
        s3 = primer1[0:i]
        s4 = primer2[-i:]
        if matched(s1, s2, i) >= continuous_matched_allow or matched(s3, s4, i) >= continuous_matched_allow:
            signal = 1
            break
    if signal == 1:
        return True # this is a cross match (reverse compliment)
    else:
        return False


def gc_tm(string):
    a = string.count("A")
    t = string.count("T")
    c = string.count("C")
    g = string.count("G")
    gc_content = (g + c) / len(string)
    tm = 2 * (a + t) + 4 * (g + c)
    return gc_content, tm


def gc_rich(string, min=0.2, max=0.8):
    gc_content, na = gc_tm(string)
    if gc_content > max or gc_content < min:
        return 1


def bad_tm(fpri, rpri, min_tm=args.min_tm, max_tm=args.max_tm,
             tm_diff=args.tm_difference, min_gc=args.min_gc, max_gc=args.max_gc):
    fgc, ftm = gc_tm(fpri)
    rgc, rtm = gc_tm(rpri)
    if fgc < min_gc or fgc > max_gc:
        return 1
    elif rgc < min_gc or rgc > max_gc:
        return 1
    elif ftm < min_tm or ftm > max_tm:
        return 1
    elif rtm < min_tm or rtm > max_tm:
        return 1
    elif abs(ftm - rtm) > tm_diff:
        return 1
    else:
        return 0

def find_target(target_pos, start, end):
    found = set()
    for i in target_pos:
        if i > start and i < end:
            found.add(i)
    return found

def make_consensus(aln):
    code = {'.': 0,
            'A': 2,
            'T': 3,
            'C': 5,
            'G': 7,
            'R': 14,
            'Y': 15,
            'M': 10,
            'K': 21,
            'S': 35,
            'W': 6 }
    code1 = {0:'.',
             2:'A',
             3:'T',
             5:'C',
             7:'G',
             14:'R',
             15:'Y',
             10:'M',
             21:'K',
             35:'S',
             6:'W'}
    snps = 0
    seq1 = aln['mat'].upper()
    seq2 = aln['pat'].upper()
    consensue = ""
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            consensue += seq1[i]
        else:
            snps += 1
            multi_base = code[seq1[i]] * code[seq2[i]]
            consensue += code1[multi_base]
    return consensue, snps


def find_extra_snp(consensue):
    pattern = re.compile(r'[RYMKSW]')
    find_all = pattern.findall(consensue)
    extra_snp = len(find_all) - 1
    return extra_snp


def get_primers(consensue, minprilen, maxprilen, stop_pos):
    pattern = re.compile(r'[RYMKSW]')
    for pos in range(stop_pos):
        for l in range(minprilen, maxprilen+1):
            tmp_pri = consensue[pos:pos+l]
            if not pattern.search(tmp_pri):
                yield tmp_pri, pos, l
            else:
                continue


def main():
    target_pos = set()
    chrome = ""
    with open(args.position, 'r') as ph:
        for i in ph.readlines():
            if i.startswith("#"):
                continue
            tmp = i.rstrip().split()
            target_pos.add(int(tmp[1]))
            chrome = tmp[0].split("_")[-1]

    fa = open(args.outpre, 'w')
    """
    log = open(args.outpre + ".primer.design.log", 'w')
    print("=", file=log)
    print("MIN_PRIMER_LEN:\t{}".format(args.min_plen), file=log)
    print("MAX_PRIMER_LEN:\t{}".format(args.max_plen), file=log)
    print("MIN_PRODUCT_LEN:\t{}".format(args.min_size), file=log)
    print("MAX_PRODUCT_LEN:\t{}".format(args.max_size), file=log)
    print("MIN_GC_CONTENT:\t{}".format(args.min_gc), file=log)
    print("MAX_GC_CONTENT:\t{}".format(args.max_gc), file=log)
    print("MIN_TM:\t{}".format(args.min_tm), file=log)
    print("MIN_TM:\t{}".format(args.max_tm), file=log)
    print("MAX_TM_DIFFERENCE:\t{}".format(args.tm_difference), file=log)
    print("MAX_CONTINUE_DIMER:\t{}".format(args.max_continous_dimer), file=log)
    print("=\n", file=log)
    """

    with open(args.infile, 'r') as fh:
        for aln in get_alignment(fh):
            alignment, phase, end = aln
            if len(target_pos) == 0:
                print("all done!")
                break
            # regard maternal genome as ref and standard
            any_target_in_this_aln = find_target(target_pos, phase['mat'], end['mat'])
            if len(any_target_in_this_aln) == 0:
                # no target snp in this alignment, so skip
                continue
            else:
                # it seems like some snp in this region, ok, let's check them one by one.
                # for example, args.max_size = 601
                for win in slide_win(alignment, phase, args.max_size):
                    tmp_aln, tmp_phase = win
                    # anchor target site at very middle position
                    wanted_pos = tmp_phase['mat'] + int((args.max_size - 1) / 2)
                    # can not allow any indel in this alignment
                    if "." in tmp_aln['mat'] or "." in tmp_aln['pat']:
                        if wanted_pos in any_target_in_this_aln:
                            #print("-\t{}\tFlanking_indel".format(wanted_pos), file=log)
                            any_target_in_this_aln.remove(wanted_pos)
                            target_pos.remove(wanted_pos)
                        continue
                    # my strategy is putting the target SNP on very middle position,
                    # then scanning upstream region and downstream to find potential 
                    # primers.
                    if wanted_pos in any_target_in_this_aln:
                        consensue, snps = make_consensus(tmp_aln)
                        consensue_rev_com = reverse_complementation(consensue)
                        scan_stop = int((args.max_size - args.min_size) / 2)
                        any_target_in_this_aln.remove(wanted_pos)
                        target_pos.remove(wanted_pos)
                        # beside the SNP of your target, there are ohter SNPs in this
                        # region, this is not what we want!
                        if find_extra_snp(consensue) > args.extra_snp_allow:
                            #print("-\t{}\t Many-snp".format(wanted_pos), file=log)
                            continue
                        else:
                            # remove this SNP which we have found from all targets, will
                            # not detect this one in next round.
                            print(">{}_{}-{}_{}\n{}".format(chrome, tmp_phase['mat'],
                                                             tmp_phase['mat']+args.max_size-2,
                                                             wanted_pos,
                                                             consensue), file=fa)
                            """
                            # find proper primer pair
                            for f in get_primers(consensue, args.min_plen,
                                                 args.max_plen, scan_stop):
                                fpri, findex, flen = f
                                real_findex = findex + 1  # 1-based
                                fend = real_findex + flen - 1

                                if ploy_x(fpri):
                                    print("F-{}:{}\tployN".format(findex, fpri), file=log)
                                    continue
                                elif gc_rich(fpri[-10:]):
                                    print("F-{}:{}\tGC_rich".format(findex, fpri), file=log)
                                    continue
                                # make a dimer with itself? also call it hairpin
                                elif dimer(fpri, fpri[::-1]):
                                    print("F-{}:{}\thairpin".format(findex, fpri), file=log)
                                    continue

                                for r in get_primers(consensue_rev_com, args.min_plen,
                                                     args.max_plen, scan_stop):
                                    # index is based on 5'->3' of reverse strand
                                    rpri, rindex, rlen = r
                                    reversed_rpri = rpri[::-1]
                                    real_rindex = args.max_size - rindex
                                    rend = real_rindex - rlen + 1
                                    product_size = real_rindex - real_findex + 1
                                    product = consensue[real_findex:real_rindex]
                                    # check poly N and 3' end, NOT GGG, NOT CCC, NOT A.
                                    if ploy_x(rpri):
                                        print("R-{}:{}\tployN".format(real_rindex, rpri), file=log)
                                        continue
                                    elif gc_rich(rpri[-10:]):
                                        print("R-{}:{}\tGC_rich".format(real_rindex, rpri), file=log)
                                        continue
                                    elif fpri[-3:] == reverse_complementation(rpri[:3]):
                                        print("F-{}:{}\tR-{}:{}\t3'_error".format(findex,
                                                                                  fpri,
                                                                                  real_rindex,
                                                                                  rpri), file=log)


                                    elif bad_tm(fpri, rpri):

                                        print("F-{}:{}\tR-{}:{}\tbam_Tm".format(findex,
                                                                                fpri,
                                                                                real_rindex,
                                                                                rpri), file=log)

                                        continue

                                    # make a dimer with itself? also call it hairpin
                                    elif dimer(rpri, rpri[::-1]):
                                        print("R-{}:{}\thairpin".format(real_rindex, rpri), file=log)
                                        continue
                                    # make a dimer with reverse primer?
                                    elif dimer(fpri, reversed_rpri):

                                        print("F-{}:{}\tR-{}:{}\tdimer".format(findex,
                                                                               fpri,
                                                                               real_rindex,
                                                                               rpri), file=log)

                                        continue
                                    else:
                                        chr_f_start = tmp_phase['mat'] + real_findex - 1
                                        chr_f_end = tmp_phase['mat'] + fend - 1
                                        chr_r_start = tmp_phase['mat'] + real_rindex - 1
                                        chr_r_end = tmp_phase['mat'] + rend -1
                                        print("{}\t{}\t{}\t{}\t{}\tF-{}\t{}-{}\t{}\t{}\t{}\tR-{}\t{}-{}\t{}".format(chrome,
                                                                                                                wanted_pos,
                                                                                                                product_size,
                                                                                                                chr_f_start,
                                                                                                                chr_f_end,
                                                                                                                real_findex,
                                                                                                                real_findex,
                                                                                                                fend, fpri,
                                                                                                                chr_r_start,
                                                                                                                chr_r_end,
                                                                                                                real_rindex,
                                                                                                                real_rindex,
                                                                                                                rend,
                                                                                                                rpri), file=fw)
                if len(any_target_in_this_aln) > 0:
                    for t in any_target_in_this_aln:
                        print("-\t{}\tEnd_Aln".format(t), file=log)
                        target_pos.remove(t)


    if len(target_pos) > 0:
        for i in target_pos:
            print("-\t{}\tNot_found".format(i), file=log)
            """
    #log.close()
    fa.close()

if __name__ == '__main__':
    main()
