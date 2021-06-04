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
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols, na=c('None'))
    bcinfo=bind_rows(bcinfo, barcodeinfo)
}

nas=bcinfo %>%
    select(everything()) %>%
    summarise_all(funs(sum(is.na(.))/dim(bcinfo)[1]))
    

contigmeans=colMeans(bcinfo[,-(1:2)])
contigmeans=aggregate(bcinfo[,-(1:2)], list(bcinfo$chrname), mean, na.rm=TRUE, na.action=NULL)

        
