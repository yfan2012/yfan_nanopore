library(tidyverse)
library(RColorBrewer)
library(cowplot)
source('barcode_plot_functions.R')

##see how filtering changes distributions, barnyard, etc

datadir='/mithril/Data/Nanopore/projects/methbin/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/barcode'

sampinfo=tibble(samp=c('neb15', 'neb17', 'neb19', 'nebdcm', 'neb11'),
                motif=c('GCNGC', 'GANTC', 'GATC', 'CCWGG', 'unmeth'))
motiflevels=seq(10, 30, 5)
threshlevels=seq(1 ,9 ,1)


####plot barcode score distributions with varying threshold filters

threshfiltpdf=file.path(dbxdir, 'threshfilt.pdf')
pdf(threshfiltpdf, h=8, w=19)
for (i in 1:dim(sampinfo)[1]) {
    threshplots=list()
    samp=sampinfo[i,]
    for (j in threshlevels) {
        bcfile=file.path(datadir, paste0(samp$samp, '_barcodes_', as.character(j), '.txt'))
        title=paste0(samp$motif, ' +.', as.character(j), ' threshold')
        plot=plot_bc_dists(bcfile, title)
        threshplots[[as.character(j)]]=plot
    }
    print(plot_grid(threshplots[[1]], threshplots[[2]], threshplots[[3]], threshplots[[4]], threshplots[[5]], threshplots[[6]], threshplots[[7]], threshplots[[8]], threshplots[[9]], ncol=3, align='v'))
}
dev.off()

summarythresh=c(1,5,9)
threshfiltsumpdf=file.path(dbxdir, 'threshfilt_summary.pdf')
pdf(threshfiltsumpdf, h=3, w=19)
for (i in 1:dim(sampinfo)[1]) {
    threshplots=list()
    samp=sampinfo[i,]
    for (j in summarythresh) {
        bcfile=file.path(datadir, paste0(samp$samp, '_barcodes_', as.character(j), '.txt'))
        title=paste0(samp$motif, ' +.', as.character(j), ' threshold')
        plot=plot_bc_dists(bcfile, title)
        threshplots[[as.character(j)]]=plot
    }
    print(plot_grid(threshplots[[1]], threshplots[[2]], threshplots[[3]], ncol=3, align='v'))
}
dev.off()




####plot motif filtered distributions

bc_cols=c('readname', 'GATC', 'GANTC', 'CCWGG', 'GCNGC')
filt_cols=c('readname', 'chrname', 'start', 'end', 'mapq')

motiffiltpdf=file.path(dbxdir, 'motiffilt.pdf')
pdf(motiffiltpdf, h=6, w=19)
for (i in 1:dim(sampinfo)[1]) {
    motifplots=list()
    samp=sampinfo[i,]
    bcfile=file.path(datadir, paste0(samp$samp, '_barcodes.txt'))
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols)
    for (j in motiflevels) {
        filtfile=file.path(datadir, paste0(samp$samp, '_barcodes_filtered_', as.character(j), '_motifs.txt'))
        filtinfo=read_tsv(filtfile, col_names=filt_cols)
        
        filtbc=barcodeinfo %>%
            filter(readname %in% filtinfo$readname)

        title=paste0(samp$motif, ' min:', as.character(j), ' motifs per read')
        plot=plot_bc_dists_fromtibble(filtbc, title)
        motifplots[[as.character(j)]]=plot
    }
    print(plot_grid(motifplots[[1]], motifplots[[2]], motifplots[[3]], motifplots[[4]], motifplots[[5]]))
}
dev.off()



summarymotifs=c(10,20,30)
motiffiltpdfsum=file.path(dbxdir, 'motiffilt_summary.pdf')
pdf(motiffiltpdfsum, h=3, w=19)
for (i in 1:dim(sampinfo)[1]) {
    motifplots=list()
    samp=sampinfo[i,]
    bcfile=file.path(datadir, paste0(samp$samp, '_barcodes.txt'))
    barcodeinfo=read_tsv(bcfile, col_names=bc_cols)
    for (j in summarymotifs) {
        filtfile=file.path(datadir, paste0(samp$samp, '_barcodes_filtered_', as.character(j), '_motifs.txt'))
        filtinfo=read_tsv(filtfile, col_names=filt_cols)
        
        filtbc=barcodeinfo %>%
            filter(readname %in% filtinfo$readname)

        title=paste0(samp$motif, ' min:', as.character(j), ' motifs per read')
        plot=plot_bc_dists_fromtibble(filtbc, title)
        motifplots[[as.character(j)]]=plot
    }
    print(plot_grid(motifplots[[1]], motifplots[[2]], motifplots[[3]], ncol=3, align='v'))
}
dev.off()


