#!/usr/bin/env python
# coding:UTF-8

from __future__ import division, print_function
from optparse import OptionParser
import sys, collections, itertools, os.path
try:
    import HTSeq
except ImportError:
    sys.stderr.write("Could not import HTSeq.")
    sys.exit(1)

def main():
    parser = OptionParser(
        usage="python %prog <in.gtf> <out.tab>",
        description="5.10 update: Add sort.\n"
                    "5.12 bug corrected: 1. No interval between exon1 and exon2. \n"
                    "2. Two genes with same gene_id but different transcript_id merge. \n"
                    "time test(hg38):\n"
                    "real    0m52.793s\n"
                    "user    0m52.170s\n"
                    "sys     0m0.596s\n",
        epilog="Author: gouge\tDate: 9/5/2017"
    )
    # parser.add_option("-g", "--gtf", action="store", dest="gtf_file", type="string",
    #                   help="input gtf annotation. ")
    (options, args) = parser.parse_args()
    gtf_file = args[0]
    out_file = args[1]
    prepare_annotation(gtf_file, out_file)

def prepare_annotation(gtf_file, out_file):
    aggregate_exons = collections.defaultdict(list)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)  # For detecting overlap
    # Load exons from GTF
    for feature in HTSeq.GFF_Reader(gtf_file, end_included=True):
    # for feature in itertools.islice(HTSeq.GFF_Reader(gtf_file, end_included=True), 500):
        if feature.type != "exon" or '_' in feature.iv.chrom:  # Filter chrom with '_'(i.e chr*alt)
            continue
        aggregate_exons[feature.attr['gene_id']].append(feature)
        exons[feature.iv] += "overlap with exon"

    # For each gene/transcript, number the exons
    for l in aggregate_exons.values():
        for i in range(len(l)):
            l[i].attr['exon_number'] = '%03d' % (i+1)

    ssList = []
    for gene in aggregate_exons.values():
        gene.sort(key=lambda f: (f.iv.chrom, f.iv.start))
        if len(gene) == 1:
            sys.stderr.write('%s has no intron. \n' % gene[0].name)
            continue
        j = 0
        for i in range(1, len(gene)):
            if gene[i-1].iv.chrom != gene[i].iv.chrom:
                continue
            if gene[i-1].iv.strand != gene[i].iv.strand:
                continue
            if gene[i-1].iv.start >= gene[i].iv.end:
                continue
            if gene[i-1].name != gene[i].name:
                continue
            if gene[i-1].attr['transcript_id'] != gene[i].attr['transcript_id']:
                continue
            j += 1
            ssList.append((gene[i-1], gene[i], j))

    # Sort and write out
    fout = open(out_file, "w")
    ssList.sort(key=lambda (g1, g2, i): (g1.iv.chrom, g1.iv.start))
    n = 1
    for gene1, gene2, i in ssList:
        fs = set()
        intron_iv = HTSeq.GenomicInterval(gene1.iv.chrom, gene1.iv.end, gene2.iv.start, gene1.iv.strand)
        for iv, s in exons[intron_iv].steps():
            fs |= s
        is_overlap = "overlap with exon" in fs
        label = ':'.join(['EIE%06d' % n, gene1.name, '%d' % i, str(is_overlap)])
        interval = gene1.iv.chrom + '_' + gene1.iv.strand + '_' +\
                   str(gene1.iv.start) + ':' + str(gene1.iv.end) + '_' + \
                   str(gene2.iv.start) + ':' + str(gene2.iv.end)
        fout.write('\t'.join([label, interval]) + '\n')
        n += 1
    fout.close()

if __name__ == '__main__':
    main()
    # prepare_annotation('hg38_ucsc.annotated.gtf', '666.tab')
