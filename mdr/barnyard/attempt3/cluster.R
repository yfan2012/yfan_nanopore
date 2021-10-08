library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')

cluster=new_cluster(7)
cluster_library(cluster, 'tidyverse')

datadir='/mithril/Data/Nanopore/projects/methbin/barnyard/barcode'
dbxdir='~/gdrive/mdr/barnyard'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)

universal=c('readname', 'chrname', 'GATC', 'GANTC', 'RAATTY', 'CATG', 'GTAC', 'SAY')

plotfile=file.path(dbxdir, '211005_mixes.pdf')
pdf(plotfile, w=13, h=7)
for (i in c(5,6)) {

    prefix=paste0('211005_mdr_barnyard_mix', as.character(i))

    barcodefile=file.path(datadir, paste0(prefix, '_barcodes.txt'))
    countsfile=file.path(datadir, paste0(prefix, '_barcodes_motifcounts.txt'))
    fullbcinfo=read_tsv(barcodefile, col_name=bc_cols, na=c('None'))

    nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
    lowna=nacount[nacount<.30]
    keepmotifs=universal
    
    bcinfo=fullbcinfo %>%
        select(all_of(keepmotifs))

    fullbccounts=read_tsv(countsfile, col_name=bc_cols, na=c('None'))
    bccounts=fullbccounts %>%
        select(all_of(keepmotifs)) %>%
        filter(across(c(-readname, -chrname), ~.x>=3))
    
    bcfilt=bcinfo %>%
        filter(readname %in% bccounts$readname)
    ecoli=bcfilt %>%
        filter(chrname=='NC_000913.3') %>%
        filter(complete.cases(.))
    staph=bcfilt %>%
        filter(chrname=='NC_007795.1') %>%
        filter(complete.cases(.))
    chrfilt=bind_rows(ecoli, staph)
    plas=bcfilt %>%
        filter(chrname=='PRW62')
    
    counts=min(c(5000, dim(ecoli)[1], dim(staph)[1]))
    bcfiltsub=bind_rows(ecoli[1:counts,], staph[1:counts,], plas)

    
    ####umap stuff
    bcdata=bcfiltsub %>%
        select(-chrname, -readname)
    
    bcumap=umap(bcdata)

    embed=tibble(x=bcumap$layout[,1],
                 y=bcumap$layout[,2],
                 label=bcfiltsub$chrname)
    embedplas=embed %>%
        filter(label=='PRW62') %>%
        mutate(shape='plas')
    embedchr=embed %>%
        filter(label!='PRW62') %>%
        mutate(shape='chr')
    
    mycolors=c(brewer.pal(8, 'Set2'), '#000000')
    myshapes=c(3,18)

    plot=ggplot(embedchr, aes(x=x, y=y, colour=label)) +
        geom_point(alpha=.5, size=.4) +
        scale_colour_manual(values=mycolors) +
        ggtitle(paste0('mix', as.character(i))) +
        theme_bw()
    mainplot=plot +
        geom_point(data=embedplas, aes(x=x, y=y, shape=shape), inherit.aes=FALSE) +
        scale_shape_manual(values=myshapes)
    print(mainplot)

    
    ####try distance stuff
    subchrfilt=bcfiltsub %>%
        filter(chrname!='PRW62')

    cluster=new_cluster(36)
    cluster_library(cluster, 'tidyverse')
    cluster_copy(cluster, 'get_read_distances')
    cluster_copy(cluster, 'classify_by_top_reads')
    cluster_copy(cluster, 'subchrfilt')
    
    plasvote=plas %>%
        rowwise() %>%
        partition(cluster)
    voteinfo=plasvote %>%
        do(classify_by_top_reads(., subchrfilt, 50))
    voteresults=voteinfo %>%
        collect() %>%
        ungroup()
    votecsvfile=file.path(datadir, paste0(prefix, '_votes.csv'))
    write_csv(voteresults, votecsvfile)
    
}
dev.off()


