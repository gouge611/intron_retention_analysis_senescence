#!/usr/bin/env python
# coding:UTF-8

from __future__ import division, print_function
from optparse import OptionParser
import sys, collections, itertools, os.path, re

try:
    import HTSeq
except ImportError:
    sys.stderr.write("Could not import HTSeq.")
    sys.exit(1)

def mapping_reads2shared_exons_introns(refGene_txt, bam_filename, minaqual, stranded, order, max_buffer_size):
    # initialise counters
    counts = {}
    counts['_empty'] = 0
    counts['_ambiguous'] = 0
    counts['_lowaqual'] = 0
    counts['_notaligned'] = 0
    counts['_ambiguous_readpair_position'] = 0

    # Read BAM file
    bam_reader = HTSeq.BAM_Reader(bam_filename)
    # CIGAR match characters (including alignment match, sequence match, and sequence mismatch
    cigar_char = ('M', '=', 'X')
    # (Refer to HTSeq-count)strand-associated
    stranded_boolean = stranded == 'yes' or stranded == 'reverse'
    reverse_boolean = stranded == 'reverse'
    def invert_strand(iv):
        iv2 = iv.copy()
        if iv2.strand == "+":
            iv2.strand = "-"
        elif iv2.strand == "-":
            iv2.strand = "+"
        else:
            raise ValueError("Illegal strand")
        return iv2

    sys.stdout.write("Gene\tfeature\trank\tposition\tlength\tread_counts\tread_counts_norm\tcoverage(%)\n")

    annot = collections.OrderedDict()
    for line in open(refGene_txt):
        gene_label, feature, rank, position, length = line.strip().split('\t')
        chrom, iv_str, strand = position.strip().split(':')
        start, end = map(int, iv_str.strip().split('-'))
        annot.setdefault(gene_label, []).append((feature, int(rank), chrom, start, end, strand, int(length)))

    for gene_name in annot:
        gene_count = {}
        gas = HTSeq.GenomicArrayOfSets("auto", stranded=stranded_boolean)
        ga = HTSeq.GenomicArray("auto", stranded=stranded_boolean, typecode="i")
        cvg_list = []

        # Annotation
        for feature, rank, chrom, start, end, strand, length in annot[gene_name]:
            iv = HTSeq.GenomicInterval(chrom, start, end, strand)
            gas[iv] += (feature, rank)
            gene_count[(feature, rank)] = 0

        # 直接对bam_reader取iter有问题，作者说是pysam的bug导致的。修正：加fetch
        boundary_left, boundary_right = min([i[3] for i in annot[gene_name]]), max([i[4] for i in annot[gene_name]])
        region_fetch = annot[gene_name][0][2] + ':' + str(int(boundary_left) - 500) + '-' + str(int(boundary_right) + 500)
        read_seq = bam_reader.fetch(region=region_fetch)

        # distinguish SE and PE mode:
        read_seq_iter = iter(bam_reader.fetch())
        one_read = next(read_seq_iter)
        pe_mode = one_read.paired_end

        if pe_mode:
            if order == 'name':
                read_seq = HTSeq.pair_SAM_alignments(read_seq)
            elif order == 'pos':
                read_seq = HTSeq.pair_SAM_alignments_with_buffer(read_seq, max_buffer_size=max_buffer_size)
            else:
                raise ValueError("Illegal order name.")

        # Mapping
        for a in read_seq:
            if not pe_mode:
                if not a.aligned:
                    counts['_notaligned'] += 1
                    continue
                if a.optional_field('NH') > 1:
                    continue
                if a.aQual < minaqual:
                    counts['_lowaqual'] += 1
                    continue
                if not reverse_boolean:
                    iv_seq = (cigop.ref_iv for cigop in a.cigar
                              if cigop.type == "M" and cigop.size > 0)
                else:
                    iv_seq = (invert_strand(cigop.ref_iv) for cigop in a.cigar
                              if cigop.type in cigar_char and cigop.size > 0)
            # pe mode
            else:
                if ((a[0] and a[0].aQual < minaqual) or
                    (a[1] and a[1].aQual < minaqual)):
                    counts['_lowaqual'] += 1
                    continue
                if ((a[0] and a[0].optional_field('NH') > 1) or
                    (a[1] and a[1].optional_field('NH') > 1)):
                    continue
                if a[0] is not None and a[0].aligned:
                    if not reverse_boolean:
                        iv_seq = (cigop.ref_iv for cigop in a[0].cigar
                                  if cigop.type in cigar_char and cigop.size > 0)
                    else:
                        iv_seq = (invert_strand(cigop.ref_iv) for cigop in a[0].cigar
                                  if cigop.type in cigar_char and cigop.size > 0)
                else:
                    iv_seq = tuple()
                if a[1] is not None and a[1].aligned:
                    if not reverse_boolean:
                        iv_seq = itertools.chain(
                            iv_seq,
                            (invert_strand(cigop.ref_iv) for cigop in a[1].cigar
                             if cigop.type in cigar_char and cigop.size > 0))
                    else:
                        iv_seq = itertools.chain(
                            iv_seq,
                            (cigop.ref_iv for cigop in a[1].cigar
                             if cigop.type in cigar_char and cigop.size > 0))

            feature_aligned = set()
            for iv in iv_seq:
                for iv2, val2 in gas[iv].steps():
                    feature_aligned |= val2
                    ga[iv] += 1  # for calculating coverage
            if len(feature_aligned) == 0:
                counts['_empty'] += 1
                continue
            # when mapping to intron, discard exons
            for f in [item for item in feature_aligned if item[0] == 'intron']:
                gene_count[f] += 1
            # when no mapping to intron, count all exons
            if 'intron' not in [x for x, y in feature_aligned]:
                for f in feature_aligned:
                    gene_count[f] += 1

        res = []
        for feature, rank, chrom, start, end, strand, length in annot[gene_name]:
            feature_count = gene_count[(feature, rank)]
            feature_count_norm = feature_count / length * 1000
            # Coverage calculation
            iv = HTSeq.GenomicInterval(chrom, start, end, strand)
            cvg_region = list(ga[iv])
            cvg = len(filter(lambda x: x > 0, cvg_region)) / len(cvg_region) * 100
            res.append([feature, rank, chrom, start, end, strand, length, feature_count, feature_count_norm, cvg])

        # Output
        for feature, rank, chrom, start, end, strand, length, feature_count, feature_count_norm, cvg in res:
            pos = "%s:%d-%d:%s" % (chrom, start, end, strand)
            sys.stdout.write('\t'.join(
                map(str, [gene_name, feature, rank, pos, length, feature_count, feature_count_norm, cvg])) + '\n')

    for fn in counts.keys():
        sys.stderr.write('%s\t%d\n' % (fn, counts[fn]))

def main():
    parser = OptionParser(
        usage="python %prog [options] <refGene.txt> <sorted.bam>",
        description='BAM file should be indexed by SAMtools. ' +
                    'Update: 1. Add coverage calculation. \n' +
                    '2. Output style adjustment. \n' +
                    '3. The way reads allocate has changed. ' +
                    'new(3.28): add paired-end mode. ',
        epilog="Author: gouge\tUpdate: 18.3.28"
    )
    parser.add_option("-a", "--minaqual", type="int", dest="minaqual",
                      default=10,
                      help="skip all reads with alignment quality lower than the given " +
                           "minimum value (default: 10)")
    parser.add_option("-s", "--stranded", type="choice", dest="stranded",
                      choices=("yes", "no", "reverse"), default="no",
                      help="'yes', 'no', or 'reverse'. Indicates whether the data is " +
                           "from a strand-specific assay (default: no ). ")
    parser.add_option("-r", "--order", dest="order",
                      choices=("pos", "name"), default="pos",
                      help="'pos' or 'name'. Sorting order of <sorted.bam> (default: pos). Paired-end sequencing " +
                            "data must be sorted either by position or by read name, and the sorting order " +
                            "must be specified. Ignored for single-end data.")
    parser.add_option("--max-reads-in-buffer", dest="max_buffer_size", type=int, default=30000000,
                      help="When <sorted.bam> is paired end sorted by position, " +
                            "allow only so many reads to stay in memory until the mates are " +
                            "found (raising this number will use more memory). Has no effect " +
                            "for single end or paired end sorted by name")
    (options, args) = parser.parse_args()

    refGene_txt = args[0]
    bam_file = args[1]

    mapping_reads2shared_exons_introns(refGene_txt,
                                       bam_file,
                                       options.minaqual,
                                       options.stranded,
                                       options.order,
                                       options.max_buffer_size)


if __name__ == '__main__':
    main()
    # mapping_reads2shared_exons_introns("annot.txt", "", 10, False)
