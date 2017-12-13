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
        description='BAM file should be indexed by SAMtools. ' +
        'Update: 1. Add coverage calculation. \n' +
        '2. Output style adjustment. \n' +
        '3. The way reads allocate has changed. ',
        epilog="Author: gouge\tUpdate: 13/Dec/2017"
    )
    parser.add_option("-a", "--minaqual", type="int", dest="minaqual",
                      default=10,
                      help="skip all reads with alignment quality lower than the given " +
                             "minimum value (default: 10)")
    parser.add_option("-s", "--stranded", type="choice", dest="stranded",
                       choices=("yes", "no", "reverse"), default="no",
                       help="'yes', 'no', or 'reverse'. Indicates whether the data is " +
                              "from a strand-specific assay (default: no ). ")
    (options, args) = parser.parse_args()

    refGene_txt = args[0]
    bam_file = args[1]
    bam_reader = HTSeq.BAM_Reader(bam_file)
    stranded = options.stranded == 'yes' or options.stranded == 'reverse'
    reverse = options.stranded == 'reverse'

    mapping_reads2shared_exons_introns(refGene_txt, bam_reader, options.minaqual, stranded)

def mapping_reads2shared_exons_introns(refGene_txt, bam_reader, minaqual=10, stranded_boolean=False):
    # initialise counters
    counts = {}
    counts[ '_empty' ] = 0
    counts[ '_ambiguous' ] = 0
    counts[ '_lowaqual' ] = 0
    counts[ '_notaligned' ] = 0
    counts['_ambiguous_readpair_position'] = 0

    sys.stdout.write("Gene\tfeature\tposition\tlength\tread_counts\tread_counts_norm\tcoverage(%)\n")

    for line in open(refGene_txt):
        gene_symbol, chrom, strand, shared_exon_start, shared_exon_end, \
        shared_intron_start, shared_intron_end, \
        sss, bbb = line.strip().split('\t')
        shared_exon_length_list, shared_intron_lengh_list = ([], [])
        gene_count = {}
        gas = HTSeq.GenomicArrayOfSets("auto", stranded = stranded_boolean)
        ga = HTSeq.GenomicArray("auto", stranded = stranded_boolean, typecode = "i")
        shared_exon_cvg_list, shared_intron_cvg_list = ([], [])

        # Annotation
        i = j = 1
        assert len(shared_exon_start.strip().split(',')) == len(shared_exon_end.strip().split(','))
        for s, e in zip(map(int, shared_exon_start.strip().split(',')), map(int, shared_exon_end.strip().split(','))):
            if s >= e:
                shared_exon_length_list.append("NA")
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            shared_exon_length_list.append(e-s)
            gas[iv] += ('exon', i)
            gene_count[('exon', i)] = 0
            i += 1
        assert len(shared_intron_start.strip().split(',')) == len(shared_intron_end.strip().split(','))
        for s, e in zip(map(int, shared_intron_start.strip().split(',')), map(int, shared_intron_end.strip().split(','))):
            if s >= e:
                shared_intron_lengh_list.append("NA")
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            shared_intron_lengh_list.append(e-s)
            gas[iv] += ('intron', j)
            gene_count[('intron', j)] = 0  # not i. T_T
            j += 1

        # Mapping
        boundary_left, boundary_right = \
            shared_exon_start.strip().split(',')[0], shared_exon_end.strip().split(',')[-1]
        for a in bam_reader.fetch(region=chrom + ':' + str(int(boundary_left) + 500) + '-' + str(int(boundary_right) + 500)):
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
                    ga[iv] += 1 # for calculating coverage
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

        # Coverage calculation
        for s, e in zip(map(int, shared_exon_start.strip().split(',')), map(int, shared_exon_end.strip().split(','))):
            if s >= e:
                shared_exon_cvg_list.append("NA")
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            cvg_region = list(ga[iv])
            cvg = len(filter(lambda x: x > 0, cvg_region)) / len(cvg_region) * 100
            shared_exon_cvg_list.append(cvg)
        for s, e in zip(map(int, shared_intron_start.strip().split(',')), map(int, shared_intron_end.strip().split(','))):
            if s >= e:
                shared_intron_cvg_list.append("NA")
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            cvg_region = list(ga[iv])
            cvg = len(filter(lambda x: x > 0, cvg_region)) / len(cvg_region) * 100
            shared_intron_cvg_list.append(cvg)

        # Output
        # sys.stdout.write(line.strip() + '\t')
        exon_count = [gene_count[fn] for fn in sorted(gene_count.keys(), key=lambda x:x[1]) if fn[0] == 'exon']
        exon_count_norm = [c/l*1000 for c, l in zip(exon_count, shared_exon_length_list)]
        intron_count = [gene_count[fn] for fn in sorted(gene_count.keys(), key=lambda x:x[1]) if fn[0] == 'intron']
        intron_count_norm = [c/l*1000 for c, l in zip(intron_count, shared_intron_lengh_list)]
        # sys.stdout.write(str(sum(intron_count)) + '\t' + ','.join(map(str, intron_count)) + '\t' + ','.join(map(str, intron_count_norm)) + '\t')
        # sys.stdout.write(str(sum(exon_count)) + '\t' + ','.join(map(str, exon_count)) + '\t' + ','.join(map(str, exon_count_norm)) + '\t')
        # sys.stdout.write(','.join(map(str, shared_intron_cvg_list)) + '\t' + ','.join(map(str, shared_exon_cvg_list)) + '\n')

        for start, end, count, count_norm, cvg in zip(map(int, shared_exon_start.strip().split(',')), map(int, shared_exon_end.strip().split(',')), exon_count, exon_count_norm, shared_exon_cvg_list):
            pos = "%s:%d-%d:%s" % (chrom, start, end, strand)
            length = end - start
            sys.stdout.write('\t'.join(map(str, [gene_symbol, "shared_exon", pos, length, count, count_norm, cvg])) + '\n')
        for start, end, count, count_norm, cvg in zip(map(int, shared_intron_start.strip().split(',')), map(int, shared_intron_end.strip().split(',')), intron_count, intron_count_norm, shared_intron_cvg_list):
            pos = "%s:%d-%d:%s" % (chrom, start, end, strand)
            length = end - start
            sys.stdout.write('\t'.join(map(str, [gene_symbol, "shared_intron", pos, length, count, count_norm, cvg])) + '\n')

    for fn in counts.keys():
        sys.stderr.write('%s\t%d\n' % (fn, counts[fn]))

if __name__ == '__main__':
    main()
    # mapping_reads2shared_exons_introns("refGene_hg19.txt_length_filter_final", HTSeq.BAM_Reader("test.bam"))
