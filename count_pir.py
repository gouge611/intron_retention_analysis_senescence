#!/usr/bin/env python
# coding:UTF-8

from __future__ import division, print_function
from optparse import OptionParser
import sys, collections, itertools, os.path, numpy, re
from scipy import stats
try:
    import HTSeq
except ImportError:
    sys.stderr.write("Could not import HTSeq.")
    sys.exit(1)

def main():
    parser = OptionParser(
        usage="python %prog [options] <EIE_junction.tab> <sorted.bam>",
        description="5.10 Basic model has finished."
                    "Note: sorted bam file is needed for input.",
#                    "\n待更：alignment异常处理、pair-end",
        epilog="Author: gouge\tDate: 10/5/2017"
    )
    parser.add_option("-a", "--minaqual", type="int", dest="minaqual",
                      default=10,
                      help="skip all reads with alignment quality lower than the given " +
                             "minimum value (default: 10)")
    parser.add_option("-s", "--stranded", type="choice", dest="stranded",
                       choices=("yes", "no", "reverse"), default="yes",
                       help="'yes', 'no', or 'reverse'. Indicates whether the data is " +
                              "from a strand-specific assay (default: yes ). " +
                              "Be sure to switch to 'no' if you use a non strand-specific RNA-Seq library " +
                              "preparation protocol. 'reverse' inverts strands and is needed for certain " +
                              "protocols, e.g. paired-end with circularization.")
    parser.add_option("-p", "--paired", type="choice", dest="paired",
                       choices=("no", "yes"), default="no",
                       help="'yes' or 'no'. Indicates whether the data is paired-end (default: no)")
    parser.add_option("-r", "--order", type="choice", dest="order",
                       choices=("pos", "name"), default="name",
                      help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
                           "data must be sorted either by position or by read name, and the sorting order " +
                           "must be specified. Ignored for single-end data." )
    (options, args) = parser.parse_args()

    ele_tab_file = args[0]
    bam_file = args[1]
    # out_file = args[2]
    stranded = options.stranded == 'yes' or options.stranded == 'reverse'
    reverse = options.stranded == 'reverse'

    junctionList = []
    for line in open(ele_tab_file):
        field = line.strip().split('\t')
        junctionList.append((field[0], field[1]))

    bam_reader = HTSeq.BAM_Reader(bam_file)
    dict_PIR = esti_each_PIR(junctionList, bam_reader, options.minaqual, stranded)

    # f = open(out_file, 'w')
    for fn in sorted(dict_PIR.keys()):
        # f.write('%s\t%d\n' % (fn, dict_PIR[fn]))
        if type(dict_PIR[fn]) == str:
            sys.stdout.write('%s\t%s\n' % (fn, dict_PIR[fn]))
        elif type(dict_PIR[fn]) == float:
            sys.stdout.write('%s\t%.2f\n' % (fn, dict_PIR[fn]))
    # f.close()

def reverse_strand(s):
    if s == '+':
        return '-'
    elif s == '-':
        return '+'
    else:
        raise  SystemError, 'illegal strand'

def test():
    14232

def esti_each_PIR(junctionELE, bam_reader, minaqual, stranded_boolean):  # 参数以元组(tuple)形式传入
    # initialize
    counts = {}
    counts['_empty'] = 0
    counts['_ambiguous'] = 0
    counts['_lowaqual'] = 0
    counts['_notaligned'] = 0
    counts['_ambiguous_readpair_position'] = 0

    pir_dict = {}
    for label, coordinates in junctionELE:
        ele_id, gene_id, num_junction, is_overlap = label.strip().split(':')
        is_overlap = is_overlap == str(True)
        chrom, strand, e1start, e1end, e2start, e2end = re.split('_|:',coordinates.strip())
        e1start, e1end, e2start, e2end = map(int, (e1start, e1end, e2start, e2end))

        features = HTSeq.GenomicArrayOfSets("auto", stranded=stranded_boolean)
        exon1_iv = HTSeq.GenomicInterval(chrom, e1start, e1end, strand)
        features[exon1_iv] += "exon1"
        exon2_iv = HTSeq.GenomicInterval(chrom, e2start, e2end, strand)
        features[exon2_iv] += "exon2"
        intron_iv = HTSeq.GenomicInterval(chrom, e1end, e2start, strand)
        features[intron_iv] += "intron"

        counts_E1 = counts_E2 = counts_E1I = counts_IE2 = counts_I = counts_E1E2 = 0
        for a in bam_reader.fetch(region=chrom + ':' + str(e1start) + '-' + str(e2end)):
            if a.iv.end > e2end:  # Discard right flanking region
                continue
            if not a.aligned:
                counts['_notaligned'] += 1
                continue
            if a.aQual < minaqual:
                counts['_lowaqual'] += 1
                continue
            if a.optional_field('NH') > 1:
                continue
            # count each mapping feature
            features_aligned = set()
            for iv, val in features[a.iv].steps():
                features_aligned |= val
            char_cigar = [cigop.type for cigop in a.cigar]
            if features_aligned == set(["exon1"]) and "N" not in char_cigar:
                counts_E1 += 1; continue
            if features_aligned == set(["exon2"]) and "N" not in char_cigar:
                counts_E2 += 1; continue
            if features_aligned == set(["intron"]) and "N" not in char_cigar:
                counts_I += 1; continue
            if features_aligned == set(["exon1", "intron"]) and "N" not in char_cigar:
                counts_E1I += 1; continue
            if features_aligned == set(["intron", "exon2"]) and "N" not in char_cigar:
                counts_IE2 += 1; continue
            if features_aligned == set(["exon1", "intron", "exon2"]) and "N" in char_cigar:
                counts_E1E2 += 1; continue

        # calculate PIR and filter
        pir = 100 * numpy.mean([counts_E1I, counts_IE2]) / (counts_E1E2 + numpy.mean([counts_E1I, counts_IE2]))
        if min(pir.item()) <= 95.0 and \
        numpy.median([counts_E1I, counts_IE2, counts_I]) + counts_E1E2 > 10.0 and \
        stats.binom_test(x=min(counts_E1I, counts_IE2, counts_I),
                         n=min(counts_E1I, counts_IE2, counts_I) + max(counts_E1I, counts_IE2, counts_I),
                         p=1/3.5, alternative="less") >= 0.05 and \
        not is_overlap:
            pir_dict[label] = pir.item()
        else:
            pir_dict[label] = 'Filtered'
    for fn in counts.keys():
        sys.stderr.write("%s\t%d\n" % (fn, counts[fn]))
    return pir_dict

if __name__ =='__main__':
    # test()
    main()
