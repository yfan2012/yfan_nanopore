library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

cluster=new_cluster(7)
cluster_library(cluster, 'tidyverse')

datadir='/mithril/Data/Nanopore/projects/methbin/barnyard/barcode'
dbxdir='~/gdrive/mdr/barnyard'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)

for (i in c(1,2,3,4)) {
    prefix=paste0('210730_mdr_barnyard_mix', as.character(i))
    barcodefile=file.path(datadir, paste0(prefix, '_barcodes.txt'))
    countsfile=file.path(datadir, paste0(prefix, '_barcodes_motifcounts.txt'))

    fullbcinfo=read_tsv(barcodefile, col_name=bc_cols, na=c('None'))

    nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
    lowna
    


