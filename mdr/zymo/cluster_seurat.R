library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'

##load info
bc_cols=c('readname', 'chrname', 'GATC', 'CCWGG', 'ATGCAT','GTCGAC','CTCCAG','CTKVAG')
barcodefile=file.path(datadir, '20190809_zymo_control_barcodes.txt')
bcinfo=read_tsv(barcodefile, col_names=bc_cols, na=c('NA', 'None')) %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    group_by(chrname) %>%
    sample_n(10000) %>%
    ungroup()
bcinfo[is.na(bcinfo)]=0

##scale
bcscaled=apply(bcinfo[,3:8], 2, scale)

##pca 
pca=prcomp(bcscaled)
pcared=pca$x[,1:4]

##umap
pcaumap=umap(pcared)
seuratfile=file.path(dbxdir, 'seurat_umap.pdf')
pdf(seuratfile, h=7, w=13)
plot=plot_umap(pcaumap$layout, bcinfo$chrname)
print(plot)
sep=plot + facet_wrap(~label)
print(sep)
dev.off()
