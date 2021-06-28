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
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols, na=c('', 'NA', 'None'))

    barcodeinfo_zeros=barcodeinfo
    barcodeinfo_zeros[is.na(barcodeinfo_zeros)]=0
    
    plot=plot_bc_dists_fromtibble(barcodeinfo_zeros, i)
    motifdists[[i]]=plot

    barcodeinfo=barcodeinfo %>%
        mutate(samp=i)
    bcinfo=bind_rows(bcinfo, barcodeinfo[1:10000,])
}

discodistspdf=file.path(dbxdir, 'dist_plots.pdf')
pdf(discodistspdf, h=6, w=19)
plot_grid(motifdists[[1]], motifdists[[2]], motifdists[[3]], motifdists[[4]], motifdists[[5]], motifdists[[6]], ncol=3, align='v')
dev.off()


####UMAP
umapplots <- function(bcinfo, bcumap, samps) {
    umaps=list()
    for (i in samps) {
        readpos=bcinfo$samp==i
        readlab=bcinfo$samp[readpos]
        readlay=bcumap$layout[readpos,]
        plot=plot_umap(readlay, readlab)
        umaps[[i]]=plot
    }
    return(umaps)
}

bcdata=bcinfo %>%
    select(-readname, -chrname, -samp)
bcdata[is.na(bcdata)]=0
bcumap=umap(bcdata)
 
umaps=umapplots(bcinfo, bcumap, samps)

umappdf=file.path(dbxdir, 'umap.pdf')
pdf(umappdf, h=10, w=13)
fullplot=plot_umap(bcumap$layout, bcinfo$samp)
print(fullplot)
subplots=fullplot +
    facet_wrap(~label)
print(subplots)
dev.off()


####pick motifs to remove manually
nacount=sort(colSums(is.na(bcinfo))/dim(bcinfo)[1])
bcreduce=bcinfo %>%
    select(-CCGCGG, -GTATAC, -GCCGGC, -GAATTC, -GCYYGAT, -ATTAAT, -RGCGCY, -GGWCC, -VGACAT, -GGTGA)
bcreduce[is.na(bcreduce)]=0

reducedata=bcreduce %>%
    select(-chrname, -samp, -readname)


##get umap and zoomed in umap
reduce_umap=umap(reducedata)
limumaps=umapplots(bcreduce, reduce_umap, samps)

reducedpdf=file.path(dbxdir, 'umap_reduced_manually.pdf')
pdf(reducedpdf, h=10, w=13)
print(plot_umap(reduce_umap$layout, bcreduce$samp))
plot_grid(limumaps[[1]], limumaps[[2]], limumaps[[3]], limumaps[[4]], limumaps[[5]], limumaps[[6]], limumaps[[7]], ncol=3, align='v')
dev.off()

##better zoomed in version
zoomdata=as_tibble(reduce_umap$layout) %>%
    mutate(samp=bcreduce$samp) %>%
    rename('x'='V1') %>%
    rename('y'='V2') 

zoompdf=file.path(dbxdir, 'umap_reduced_manuallylim.pdf')
pdf(zoompdf, h=10, w=13)
zoomplot=plot_umap(reduce_umap$layout, bcreduce$samp, c(-10, 10), c(-10, 10))
print(zoomplot)
zoomplot_sep=zoomplot +
    facet_wrap(~label)
print(zoomplot_sep)
dev.off()


bcreduce2=bcreduce %>%
    select(-CCWGG, -GTAC, -TCGA, -CCGG)
reducedata2=bcreduce2 %>%
    select(-chrname, -samp, -readname)

reduce_umap2=umap(reducedata2)

zoom2pdf=file.path(dbxdir, 'umap_reduced_manuallylim2.pdf')
pdf(zoom2pdf, h=10, w=13)
zoomplot=plot_umap(reduce_umap2$layout, bcreduce2$samp, c(-15, 15), c(0,25))
print(zoomplot)
zoomplot_sep=zoomplot +
    facet_wrap(~label)
print(zoomplot_sep)
dev.off()


