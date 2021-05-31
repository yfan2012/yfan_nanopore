library(tidyverse)
library(cowplot)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/disco/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/disco'

samps=list.dirs(datadir, full.names=FALSE, recursive=FALSE)

barcodesfile='~/Code/yfan_nanopore/mdr/disco/disco_barcodes.txt'
barcodes=read_tsv(barcodesfile, col_names=c('motif'))

bc_cols=c('readname', 'chrname', barcodes$motif)

bcinfo=NULL
motifdists=list()

for (i in samps) {
    bcfile=file.path(datadir, i, paste0(i, '_barcodes.txt'))
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols)

    plot=plot_bc_dists_fromtibble(barcodeinfo, i)
    motifdists[[i]]=plot

    barcodeinfo=barcodeinfo %>%
        mutate(samp=i)
    bcinfo=bind_rows(bcinfo, barcodeinfo[1:10000,])
}

discodistspdf=file.path(dbxdir, 'dist_plots.pdf')
pdf(discodistspdf, h=6, w=19)
plot_grid(motifdists[[1]], motifdists[[2]], motifdists[[3]], motifdists[[4]], motifdists[[5]], motifdists[[6]], ncol=3, align='v')
dev.off()


bcdata=bcinfo %>%
    select(-readname, -chrname, -samp)

bcumap=umap(bcdata)
umaps=list()

for (i in samps) {
    readpos=bcinfo$samp==i
    readlab=bcinfo$samp[readpos]
    readlay=bcumap$layout[readpos,]
    plot=plot_umap(readlay, readlab)
    umaps[[i]]=plot
}

umappdf=file.path(dbxdir, 'umap.pdf')
pdf(umappdf, h=7, w=13)
print(plot_umap(bcumap$layout, bcinfo$samp))
plot_grid(umaps[[1]], umaps[[2]], umaps[[3]], umaps[[4]], umaps[[5]], umaps[[6]], ncol=3, align='v')
dev.off()
    
