library(tidyverse)
library(cowplot)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/disco/barcode'
dbxdir='~/gdrive/mdr/disco'
i='MinION_JM3O_NAT'

##read in barcode data
barcodesfile='~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt'
barcodes=read_tsv(barcodesfile, col_names=c('motif'))
bc_cols=c('readname', 'chrname', barcodes$motif)

bcfile=file.path(datadir, i, paste0(i, '_barcodes.txt'))
bcinfo=read_tsv(bcfile, col_names=bc_cols, na=c('', 'NA', 'None'))

##throw out certain motifs and incomplete reads
nas=colSums(is.na(bcinfo))/dim(bcinfo)[1]
keepmotifs=names(nas[nas<.3])

bcfilt=bcinfo %>%
    select(all_of(keepmotifs)) %>%
    filter(complete.cases(.))
bcdata=bcfilt %>%
    select(-readname, -chrname)


##subsamp to 5k reads max per ref genome
chrcounts=table(bcfilt$chrname)
keepchrs=names(-sort(-chrcounts)[1:12])
bcfiltsub=bcfilt %>%
    rowwise() %>%
    filter(chrname %in% keepchrs) %>%
    group_by(chrname) %>%
    sample_n(1000) %>%
    ungroup()
bcdata=bcfiltsub %>%
    select(-readname, -chrname)

##umap
bcumap=umap(bcdata)

##plot
plotfile=file.path(dbxdir, 'mouse_ref_umap.pdf')
pdf(plotfile, h=7, w=13)
plot=plot_umap(bcumap$layout, bcfiltsub$chrname)
print(plot)
dev.off()



