library(tidyverse)
library(cowplot)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/disco/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/disco'

samps=c('BA', 'BF', 'CP', 'HP', 'MH', 'NG', 'TP')

barcodesfile='~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt'
barcodes=read_tsv(barcodesfile, col_names=c('motif'))

bc_cols=c('readname', 'chrname', barcodes$motif)

bcinfo=NULL
motifdists=list()


for (sp in samps) {
    i=paste0('MinION_', sp, '_NAT')
    bcfile=file.path(datadir, i, paste0(i, '_barcodes.txt'))
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols, na=c('', 'NA', 'None'))

    barcodeinfo=barcodeinfo %>%
        mutate(samp=i)
    bcinfo=bind_rows(bcinfo, barcodeinfo[1:30000,])
}

numreads=min(table(bcinfo$chrname))

nacounts=colSums(is.na(bcinfo))/dim(bcinfo)[1]
