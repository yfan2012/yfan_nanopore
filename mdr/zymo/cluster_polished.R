library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/mdr/cluster_functions.R')


datadir='/mithril/Data/Nanopore/projects/methbin/zymo'
dbxdir='~/gdrive/mdr/mdr'
prefix='20190809_zymo_control_polished'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)


barcodedatafile=file.path(datadir, 'barcode', prefix, '20190809_zymo_control_contigs_barcodes15.txt')
fullbcinfo=read_tsv(barcodedatafile, col_names=bc_cols, na=c('None'))
countsfile=file.path(datadir, 'barcode', prefix, '20190809_zymo_control_barcodes15_motifcounts.txt')
fullbccounts=read_tsv(countsfile, col_name=bc_cols, na=c('None'))


nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
lowna=nacount[nacount<.2]
keepmotifs=names(lowna)

bcinfo=fullbcinfo %>%
    select(all_of(keepmotifs))
bccounts=fullbccounts %>%
    select(all_of(keepmotifs)) %>%
    filter(across(c(-readname, -chrname), ~.x>=5))

bcfilt=bcinfo %>%
    filter(readname %in% bccounts$readname)


