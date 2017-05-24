#!/usr/bin/env python
# coding:UTF-8

from __future__ import division, print_function
from optparse import OptionParser
import sys, collections, itertools, os.path, numpy, re
try:
    import HTSeq
except ImportError:
    sys.stderr.write("Could not import HTSeq.")
    sys.exit(1)


def main():
    parser = OptionParser(
        usage="python %prog [options] <refGene.txt> <sorted.bam>",
        description='BAM file should be indexed. ',
        epilog="Author: gouge\tDate: 18/5/2017\tFIxed on 24/5/2017"
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
    (options, args) = parser.parse_args()

    refGene_txt = args[0]
    bam_file = args[1]
    bam_reader = HTSeq.BAM_Reader(bam_file)
    stranded = options.stranded == 'yes' or options.stranded == 'reverse'
    reverse = options.stranded == 'reverse'

    mapping_reads2shared_exons_introns(refGene_txt, bam_reader, options.minaqual, stranded)

def mapping_reads2shared_exons_introns(refGene_txt, bam_reader, minaqual=10, stranded_boolean=True):
    # initialise counters
    counts = {}
    counts[ '_empty' ] = 0
    counts[ '_ambiguous' ] = 0
    counts[ '_lowaqual' ] = 0
    counts[ '_notaligned' ] = 0
    counts['_ambiguous_readpair_position'] = 0

    for line in open(refGene_txt):
        gene_symbol, chrom, strand, shared_exon_start, shared_exon_end, \
        shared_intron_start, shared_intron_end, \
        sss, bbb = line.strip().split('\t')
        shared_exon_length_list = shared_intron_lengh_list = []
        gene_count = {}
        gas = HTSeq.GenomicArrayOfSets("auto", stranded=stranded_boolean)

        i = j = 1
        assert len(shared_exon_start.strip().split(',')) == len(shared_exon_end.strip().split(','))
        for s, e in zip(map(int, shared_exon_start.strip().split(',')), map(int, shared_exon_end.strip().split(','))):
            if s >= e:
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            shared_exon_length_list.append(e-s)
            gas[iv] += ('exon', i)
            gene_count[('exon', i)] = 0
            i += 1
        assert len(shared_intron_start.strip().split(',')) == len(shared_intron_end.strip().split(','))
        for s, e in zip(map(int, shared_intron_start.strip().split(',')), map(int, shared_intron_end.strip().split(','))):
            if s >= e:
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            shared_intron_lengh_list.append(e-s)
            gas[iv] += ('intron', j)
            gene_count[('intron', j)] = 0  # not i. T_T
            j += 1

        boundary_left, boundary_right = shared_exon_start.strip().split(',')[0], shared_exon_end.strip().split(',')[-1]
        for a in bam_reader.fetch(region=chrom + ':' + str(boundary_left) + '-' + str(boundary_right)):
            if not a.aligned:
                counts['_notaligned'] += 1
                continue
            if a.optional_field('NH') > 1:
                continue
            if a.aQual < minaqual:
                counts['_lowaqual'] += 1
                continue

            feature_aligned = set()
            for cigop in a.cigar:
                if cigop.type != 'M':
                    continue
                for iv, val in gas[cigop.ref_iv].steps():
                    feature_aligned |= val
            if len(feature_aligned) == 0:
                counts['_empty'] += 1
                continue
            for f in [item for item in feature_aligned if item[0] == 'intron']:
                # when mapping to intron, discard exons
                gene_count[f] += 1
            if 'intron' not in [x for x, y in feature_aligned]:
                # when no mappint to intron, count all exons
                for f in feature_aligned:
                    gene_count[f] += 1

        # Output
        sys.stdout.write(line.strip() + '\t')
        exon_count = [gene_count[fn] for fn in sorted(gene_count.keys(), key=lambda x:x[1]) if fn[0] == 'exon']
        exon_count_norm = [c/l for c, l in zip(exon_count, shared_exon_length_list)]
        intron_count = [gene_count[fn] for fn in sorted(gene_count.keys(), key=lambda x:x[1]) if fn[0] == 'intron']
        intron_count_norm = [c/l for c, l in zip(intron_count, shared_intron_lengh_list)]
        sys.stdout.write(str(sum(intron_count)) + '\t' + ','.join(map(str, intron_count)) + '\t' + ','.join(map(str, intron_count_norm)) + '\t')
        sys.stdout.write(str(sum(exon_count)) + '\t' + ','.join(map(str, exon_count)) + '\t' + ','.join(map(str, exon_count_norm)) + '\n')
    for fn in counts.keys():
        sys.stderr.write('%s\t%d\n' % (fn, counts[fn]))

if __name__ == '__main__':
    main()
    # mapping_reads2shared_exons_introns("refGene_hg19.txt_length_filter_final", HTSeq.BAM_Reader("test.bam"))
