library(tidyverse)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/disco/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/disco'

samps=list.dirs(datadir, full.names=FALSE, recursive=FALSE)

barcodesfile='~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt'
barcodes=read_tsv(barcodesfile, col_names=c('motif'))
bc_cols=c('readname', 'chrname', barcodes$motif)

bcinfo=NULL
for (i in samps) {
    bcfile=file.path(datadir, i, paste0(i, '_barcodes.txt'))
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols)
    bcinfo=bind_rows(bcinfo, barcodeinfo)
}
    
contigmeans=aggregate(bcinfo[,-(1:2)], list(bcinfo$chrname), mean)

        
