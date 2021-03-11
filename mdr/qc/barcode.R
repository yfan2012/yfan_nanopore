library(tidyverse)
library(RColorBrewer)

datadir='~/data/mdr/qc/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/barcode'
sampinfo=tibble(samp=c('neb15', 'neb17', 'neb19', 'nebdcm', 'neb11'),
                motif=c('GCNGC', 'GANTC', 'GATC', 'CCWGG', 'unmeth'))


### for plotting dists
bc_cols=c('readname', 'GATC', 'GANTC', 'CCWGG', 'GCNGC')
barcodeinfo=tibble(readname=as.character(),
                   samp=as.character(),
                   motif=as.character(),
                   count=as.numeric())

bcdistsfile=file.path(dbxdir, 'score_dists.pdf')
pdf(bcdistfile, h=8, w=15)
for (i in 1:dim(sampinfo)[1]) {
    samp=sampinfo[i,]
    bcfile=file.path(datadir, paste0(samp$samp, '_barcodes.txt'))
    bc=read_tsv(bcfile, col_names=bc_cols) %>%
        mutate(samp=samp$motif) %>%
        gather('motif', 'count', -readname, -samp)
    barcodeinfo=bind_rows(barcodeinfo, bc)
    plot=ggplot(bc, aes(x=count, colour=motif, fill=motif, alpha=.2)) +
        geom_density() +
        ggtitle(paste0(samp$samp, '_', samp$motif)) +
        xlab('Normalized Counts') +
        scale_x_log10() +
        scale_fill_brewer(palette = 'Set2') +
        scale_colour_brewer(palette = 'Set2') +
        theme_bw()
    print(plot)
}
dev.off()
    

### for pca
library(ggbiplot)
bcs=tibble(readname=as.character(),
           GATC=as.numeric(),
           GANTC=as.numeric(),
           CCWGG=as.numeric(),
           GCNGC=as.numeric(),
           samp=as.character())
for (i in 1:dim(sampinfo)[1]) {
    samp=sampinfo[i,]
    bcfile=file.path(datadir, paste0(samp$samp, '_barcodes.txt'))
    bc=read_tsv(bcfile, col_names=bc_cols) %>%
        mutate(samp=samp$motif)
    bcs=bind_rows(bcs, bc[1:500,])
}
pca=prcomp(bcs[,c(2:5)], center=TRUE, scale.=TRUE)
pcapdf=file.path(dbxdir, 'pca_test.pdf')
pdf(pcapdf, h=10, w=6)
print(ggbiplot(pca, groups=bcs$samp) +
      scale_color_brewer(palette = 'Set2') +
      theme_bw() +
      theme(legend.direction = 'horizontal', legend.position = 'top')) 
dev.off()

