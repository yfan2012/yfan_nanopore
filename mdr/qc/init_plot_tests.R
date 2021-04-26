library(tidyverse)
library(RColorBrewer)
library(cowplot)
source('barcode_plot_functions.R')

##test some plotting functions initially


##datadir='~/data/mdr/qc/barcode'
datadir='/mithril/Data/Nanopore/projects/methbin/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/barcode'
sampinfo=tibble(samp=c('neb15', 'neb17', 'neb19', 'nebdcm', 'neb11'),
                motif=c('GCNGC', 'GANTC', 'GATC', 'CCWGG', 'unmeth'))

### for plotting barnyard
bcfile1=file.path(datadir, paste0(sampinfo[2,]$samp, '_barcodes.txt'))
bcfile2=file.path(datadir, paste0(sampinfo[3,]$samp, '_barcodes.txt'))
abarnplot=plot_barnyard(bcfile1, bcfile2, 'GATC', 'GANTC')

cbcfile1=file.path(datadir, paste0(sampinfo[1,]$samp, '_barcodes.txt'))
cbcfile2=file.path(datadir, paste0(sampinfo[4,]$samp, '_barcodes.txt'))
cbarnplot=plot_barnyard(cbcfile1, cbcfile2, 'GCNGC', 'CCWGG')

barnyardfile=file.path(dbxdir, 'barnyard.pdf')
pdf(barnyardfile, h=8, w=11)
print(abarnplot)
print(cbarnplot)
dev.off()



### for plotting dists
bc_cols=c('readname', 'GATC', 'GANTC', 'CCWGG', 'GCNGC')
bcplots=list()
singledistsfile=file.path(dbxdir, 'score_dists_single.pdf')
pdf(singledistsfile, h=8, w=15)
for (i in 1:dim(sampinfo)[1]) {
    samp=sampinfo[i,]
    bcfile=file.path(datadir, paste0(samp$samp, '_barcodes.txt'))
    plot=plot_bc_dists(samp)
    print(plot)
    bcplots[[samp$samp]]=plot
}
dev.off()

bcdistfile=file.path(dbxdir, 'score_dists.pdf')
pdf(bcdistfile, h=10, w=18)
plot_grid(bcplots[[2]], bcplots[[1]], bcplots[[3]], bcplots[[4]], bcplots[[5]], ncol=2, align='v')
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



### for zymo
zymobcs=tibble(readname=as.character(),
               GATC=as.numeric(),
               GANTC=as.numeric(),
               CCWGG=as.numeric(),
               GCNGC=as.numeric(),
               samp=as.character())
zymobcfile=file.path('~/data/mdr/zymo/barcode/zymo_barcodes_sub.txt')
zymobc=read_tsv(zymobcfile, col_names=bc_cols)
zymopca=prcomp(zymobc[,c(2:5)], center=TRUE, scale.=TRUE)
zymopcapdf=file.path(dbxdir, 'pca_zymo.pdf')
pdf(zymopcapdf)
print(ggbiplot(zymopca) +
      theme_bw() +
      theme(legend.direction = 'horizontal', legend.position = 'top'))
dev.off()


