# Intron Retention Analysis of RNA-seq data during cellular senescence

This repository contains in-house scripts that are used for intron retention analysis in the paper
- "Prevalent intron retention fine-tunes gene expression and contributes to cellular senescence" by Jun Yao, Dong Ding, Xueping Li, Ting Shen, Haihui Fu, Hua Zhong, Gang Wei, Ting Ni. 

## Description of files

`get_IR_annotation.py` is for obtaining potential IR events annotation from refGene downloaded from UCSC Genome Browser database. Annotation file in hg19 can be downloaded from `ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz`. You can also change `hg19` in the url to any other assembly like `hg38` or `mm9`. 

`shared_exon_intron_mapped_reads.py`is to generate a count file contained read count in intron and flanking exon(s) for each potential IR event. 

`calc_IRI.R` is to calculate IRI (Intron Retention Index) for each potential IR event based on the read count file generated in last step. 
