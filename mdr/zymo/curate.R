library(tidyverse)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')

##curated barcodes
projdir='/mithril/Data/Nanopore/projects/methbin/zymo'
datadir=file.path(projdir, 'curate')
prefix='20190809_zymo_control'
dbxdir='~/gdrive/mdr/zymo'


##read in barcodes
motiffile='~/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt'
motifinfo=read_tsv(motiffile, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)


##read in data
barcodefile=file.path(datadir, paste0(prefix, '_curated.txt'))
bcinfo=read_tsv(barcodefile, col_names=bc_cols, na=c('NA', 'None'))
countsfile=file.path(datadir, paste0(prefix, '_motifcounts_curated.txt'))
bccounts=read_tsv(countsfile, col_names=bc_cols, na=c('NA', 'None'))


nacount=colSums(is.na(bcinfo)/dim(bcinfo)[1])


countsfilter=bccounts %>%
    select(-AGCCGCC, -CTCGAG, -GTCGAC, -TCTAGA, -CCCGGG, -GGTCTC, -CTGGAG, -TGGCCA, -CTGCAG) %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    filter(across(c(-readname, -chrname), ~ .x>=2))
