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
        description="5.10 update:\n"
                    "Add sort.\n"
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
    # Load exons from GTF
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    for feature in HTSeq.GFF_Reader(gtf_file, end_included=True):
    # for feature in itertools.islice(HTSeq.GFF_Reader(gtf_file, end_included=True), 500):
        if feature.type == "exon":
            exons[feature.iv] += (feature.attr['gene_id'], feature.attr['transcript_id'])

    # Aggregate all exons to each gene/transcript
    aggregate_exons = collections.defaultdict(list)
    for iv, s in exons.steps():
        if len(s) == 0:
            continue
        gene_id = list(s)[0][0]
        f = HTSeq.GenomicFeature(gene_id, "exon", iv)
        f.source = os.path.basename(sys.argv[0])
        f.attr = {}
        f.attr['gene_id'] = gene_id
        f.attr['transcript_id'] = list(s)[0][1]
        aggregate_exons[gene_id].append(f)
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
        for i in range(1, len(gene)):
            if gene[i-1].iv.chrom != gene[i].iv.chrom:
                continue
            if gene[i-1].iv.strand != gene[i].iv.strand:
                continue
            if gene[i-1].iv.start >= gene[i].iv.end:
                continue
            if gene[i-1].name != gene[i].name:
                continue
            ssList.append((gene[i-1], gene[i], i))

    # Sort and write out
    fout = open(out_file, "w")
    ssList.sort(key=lambda (g1, g2, i): (g1.iv.chrom, g1.iv.start))
    n = 1
    for gene1, gene2, i in ssList:
        label = ':'.join(['EIE%06d' % n, gene1.name, '%d' % i])
        interval = gene1.iv.chrom + '_' + gene1.iv.strand + '_' +\
                   str(gene1.iv.start) + ':' + str(gene1.iv.end) + '_' + \
                   str(gene2.iv.start) + ':' + str(gene2.iv.end)
        fout.write('\t'.join([label, interval]) + '\n')
        n += 1
    fout.close()

if __name__ == '__main__':
    # main()
    prepare_annotation('hg38_ucsc.annotated.gtf', '666.gtf')
