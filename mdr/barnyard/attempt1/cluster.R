library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

##realized that the plasmids in these mixes were lost in staph
##plasmids in ecoli were retained, but lost a huge insert
##just check to see how well the chr reads separate



cluster=new_cluster(7)
cluster_library(cluster, 'tidyverse')

datadir='/mithril/Data/Nanopore/projects/methbin/barnyard/barcode'
dbxdir='~/gdrive/mdr/barnyard'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)



plotfile=file.path(dbxdir, '210730_mixes.pdf')
pdf(plotfile, w=13, h=7)
for (i in c(1,2,3,4)) {
    prefix=paste0('210730_mdr_barnyard_mix', as.character(i))
    barcodefile=file.path(datadir, paste0(prefix, '_barcodes.txt'))
    countsfile=file.path(datadir, paste0(prefix, '_barcodes_motifcounts.txt'))

    fullbcinfo=read_tsv(barcodefile, col_name=bc_cols, na=c('None'))

    nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
    lowna=nacount[nacount<.3]
    keepmotifs=names(lowna)

    bcinfo=fullbcinfo %>%
        select(all_of(keepmotifs))

    fullbccounts=read_tsv(countsfile, col_name=bc_cols, na=c('None'))
    bccounts=fullbccounts %>%
        select(all_of(keepmotifs)) %>%
        filter(across(c(-readname, -chrname), ~.x>=5))
    
    bcfilt=bcinfo %>%
        filter(readname %in% bccounts$readname)
    
    ecoli=bcfilt %>%
        filter(chrname=='NC_000913.3') %>%
        filter(complete.cases(.))
    staph=bcfilt %>%
        filter(chrname=='NC_007795.1') %>%
        filter(complete.cases(.))

    counts=min(c(5000, dim(ecoli)[1], dim(staph)[1]))
    
    bcfiltsub=bind_rows(ecoli[1:counts,], staph[1:counts,])

    bcdata=bcfiltsub %>%
        select(-chrname, -readname)

    bcumap=umap(bcdata)

    plot=plot_umap(bcumap$layout, bcfiltsub$chrname, alpha=.7) +
        ggtitle(paste0('mix', as.character(i)))
    print(plot)
}
dev.off()


