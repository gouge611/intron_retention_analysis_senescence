from __future__ import division, print_function

import collections
import copy
import sys
import itertools
from optparse import OptionParser

try:
    import HTSeq
except ImportError:
    sys.stderr.write("Could not import HTSeq.")
    sys.exit(1)

def main():
    parser = OptionParser(
        usage="python %prog <refGene.txt> <output.prefix>",
        description= "Deal with RefGene model(genePred) downloaded from UCSC. "
                     "Output: <prefix>_annot.txt describes shared exons and introns; "
                     "<prefix>_join.txt describes introns neibored with exon(s). ",
        epilog="Author: gouge"
    )
    (options, args) = parser.parse_args()

    refgene = args[0]
    output_prefix = args[1]

    Anno = genesymbol_convert(refgene)
    # get_overlap_region(Anno, output_prefix, chromsizes)
    get_overlap_region_new(Anno, output_prefix)

def genesymbol_convert(refGene_filename):
    gene_label_dict = collections.OrderedDict()
    for line in open(refGene_filename):
    # for line in itertools.islice(open(refGene_filename), 366, 369): # test overlap gene pair: CENPS/CENPS-CORT
        # colnames() = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames") # in R
        bin, name, chrom, strand, \
        txStart, txEnd, cdsStart, cdsEnd, \
        exonCount, exonStarts, exonEnds, \
        score, name2,\
        cdsStartStat, cdsEndStat, exonFrames = line.strip().split('\t')
        if chrom.find('_') != -1:
            continue

        label = "/".join([chrom, name2, strand])
        fields = [bin, name, chrom, strand, int(txStart), int(txEnd), cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames]
        if label in gene_label_dict:
            gene_label_dict[label].append(fields)
        else:
            gene_label_dict[label] = [fields]

    def getOverlapFeatures(key, val):
        assert isinstance(val, list)
        ls_sorted = sorted(val, key=lambda x:(x[4], x[5]))
        n = 1
        ls_sorted[0].append(key + "/" + str(n))
        for i in range(1, len(ls_sorted)):
            if ls_sorted[i-1][5] > ls_sorted[i][4]:
                ls_sorted[i].append(key + "/" + str(n))
            else:
                n += 1
                ls_sorted[i].append(key + "/" + str(n))
        return ls_sorted

    converted_list = []
    for genename in gene_label_dict:
        gene_label_dict[genename] = getOverlapFeatures(genename, gene_label_dict[genename])
        for ls in gene_label_dict[genename]:
            converted_list.append(tuple(ls))
    return converted_list

def get_overlap_region(anno, filename_output_prefix, filename_chromsizes):
    gas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    for line in anno:
        bin, name, chrom, strand, \
        txStart, txEnd, cdsStart, cdsEnd, \
        exonCount, exonStarts, exonEnds, \
        score, name2,\
        cdsStartStat, cdsEndStat, exonFrames, label = line

        exonStartsList = map(int, exonStarts.strip().rstrip(',').split(','))
        exonEndsList = map(int, exonEnds.strip().rstrip(',').split(','))
        intronStartList = copy.deepcopy(exonEndsList)[1:len(exonEndsList)-1]
        intronEndList = copy.deepcopy(exonStartsList)[2:len(exonStartsList)]

        for s, e in zip(exonStartsList, exonEndsList):
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            gas[iv] += (label, chrom, "exon", strand)
        for s, e in zip(intronStartList, intronEndList):
            if e - s < 1:
                continue
            iv = HTSeq.GenomicInterval(chrom, s, e, strand)
            gas[iv] += (label, chrom, "intron", strand)
    # return gas

    tmp_label = ''
    gene_dict = collections.OrderedDict()
    file_annot = open(filename_output_prefix + "_annot.txt", "w")
    for chr, chr_len in [line.strip().split('\t') for line in open(filename_chromsizes)]:
        for chrIV in [HTSeq.GenomicInterval(chr, 0, int(chr_len), strand) for strand in ("+", "-")]:
            for iv, val in gas[chrIV].steps():
                # if len(val) <= 1:
                #     continue
                # print(iv, val, file=file_annot)
                if len(val) != 1:
                    continue
                if iv.end - iv.start < 10:
                    continue
                if tuple(val)[0][0] == tmp_label:
                    n += 1
                else:
                    n = 1
                    tmp_label = tuple(val)[0][0]
                pos = ('%s:%d-%d:%s') % (iv.chrom, iv.start, iv.end, iv.strand)
                gene_dict.setdefault(tuple(val)[0][0], []).append((tuple(val)[0][0], tuple(val)[0][2], n, pos, iv.end - iv.start))
                # print('\t'.join(map(str, [tuple(val)[0][0], tuple(val)[0][2], n, pos, iv.end - iv.start])))
                file_annot.write('\t'.join(map(str, [tuple(val)[0][0], tuple(val)[0][2], n, pos, iv.end - iv.start])) + '\n')
    file_annot.close()

    file_join = open(filename_output_prefix + '_join.txt', 'w')
    for gene_name in gene_dict:
        ls_exon = [x[2] for x in gene_dict[gene_name] if x[1] == "exon"]
        ls_intron = [x for x in gene_dict[gene_name] if x[1] == "intron"]
        for y in ls_intron:
            try:
                exon_rank_left = max([x for x in ls_exon if x < y[2]])
            except:
                exon_rank_left = -1
            try:
                exon_rank_right = min([x for x in ls_exon if x > y[2]])
            except:
                exon_rank_right = -1
            file_join.write('\t'.join(map(str, y)) + '\t' + str(exon_rank_left) + '\t' + str(exon_rank_right) + '\n')
    file_join.close()
def get_overlap_region_new(anno, filename_output_prefix):
    # allocate transcripts to each gene
    gene_dict = collections.OrderedDict()
    for line in anno:
        bin, name, chrom, strand, \
        txStart, txEnd, cdsStart, cdsEnd, \
        exonCount, exonStarts, exonEnds, \
        score, name2,\
        cdsStartStat, cdsEndStat, exonFrames, label = line
        gene_dict.setdefault(label, []).append(line)

    # put all transcripts into local gas, and get shared features for each gene
    shared_feature_dict = collections.OrderedDict()
    # and then put all features in global gas
    gas_global = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    
    for label in gene_dict:
        gas_local = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        for line in gene_dict[label]:
            bin, name, chrom, strand, \
            txStart, txEnd, cdsStart, cdsEnd, \
            exonCount, exonStarts, exonEnds, \
            score, name2,\
            cdsStartStat, cdsEndStat, exonFrames, label = line
            exonStartsList = map(int, exonStarts.strip().rstrip(',').split(','))
            exonEndsList = map(int, exonEnds.strip().rstrip(',').split(','))
            intronStartList = copy.deepcopy(exonEndsList)[0:len(exonEndsList)-1]
            intronEndList = copy.deepcopy(exonStartsList)[1:len(exonStartsList)]
            for s, e in zip(exonStartsList, exonEndsList):
                iv = HTSeq.GenomicInterval(chrom, s, e, strand)
                gas_local[iv] += (label, "exon")
                gas_global[iv] += (label, "exon")
            for s, e in zip(intronStartList, intronEndList):
                if e - s < 1:
                    continue
                iv = HTSeq.GenomicInterval(chrom, s, e, strand)
                gas_local[iv] += (label, "intron")
                gas_global[iv] += (label, "intron")
        boundary_left, boundary_right = min([x[4] for x in gene_dict[label]]), max([x[5] for x in gene_dict[label]])
        for iv, val in gas_local[HTSeq.GenomicInterval(chrom, boundary_left, boundary_right, strand)].steps():
            if len(val) == 1:
                shared_feature_dict.setdefault(tuple(val)[0][0], []).append((iv, tuple(val)[0][1]))

    # check whether a feature overlaps with other genes, and return the state
    def check_overlap(iv, gas):
        assert isinstance(iv, HTSeq.GenomicInterval)
        assert isinstance(gas, HTSeq.GenomicArrayOfSets)
        features = set()
        for iv2, val in gas[iv].steps():
            features |= val
        if len(features) == 1:
            return 'Clean'
        elif len(features) > 1:
            return 'Overlap'
        else:
            return 'NA'

    join_dict = collections.OrderedDict()
    file_annot = open(filename_output_prefix + "_annot.txt", "w")
    file_join = open(filename_output_prefix + '_join.txt', 'w')
    for label in shared_feature_dict:
        # filter features overlaping with other genes
        features_filtered = [item for item in shared_feature_dict[label] if check_overlap(item[0], gas_global) == 'Clean']
        for iv, feature in features_filtered:
            rank = features_filtered.index((iv, feature)) + 1
            pos = ('%s:%d-%d:%s') % (iv.chrom, iv.start, iv.end, iv.strand)
            # join_dict.setdefault(label, []).append((label, feature, rank, pos, iv.end - iv.start))
            file_annot.write('\t'.join(map(str, [label, feature, rank, pos, iv.end - iv.start])) + '\n')

            if feature == 'intron':
                rank_exons = [features_filtered.index(x) + 1 for x in features_filtered if x[1] == "exon"]
                try:
                    exon_rank_left = max([x for x in rank_exons if x < rank])
                except:
                    exon_rank_left = -1
                try:
                    exon_rank_right = min([x for x in rank_exons if x > rank])
                except:
                    exon_rank_right = -1
                file_join.write('\t'.join(map(str, [label, feature, rank, pos, iv.end - iv.start, exon_rank_left, exon_rank_right])) + '\n')
    file_annot.close()
    file_join.close()

if __name__ == '__main__':
    import time
    start = time.clock()
    main()
    # Anno = genesymbol_convert('refGene.txt')
    # get_overlap_region_new(Anno, 'hg19')
    end = time.clock()
    print(end-start)