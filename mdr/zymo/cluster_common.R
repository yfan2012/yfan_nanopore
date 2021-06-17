library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

cluster=new_cluster(7)
cluster_library(cluster, 'tidyverse')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/barcode'
srcdir='~/Code/yfan_nanopore/mdr/rebase'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'

conds=c('10', '15', '20')

checkfilt <- function(chrblock, num, countsfilter) {
    ##subsample reads from chr
    blockfilt=countsfilter %>%
        filter(chrname==chrblock$chrname[1])
    sub=tibble()
    read=1
    while (dim(sub)[1]<num && read<dim(chrblock)[1]) {
        readinfo=chrblock[read,]
        read=read+1
        if (readinfo$readname %in% blockfilt$readname) {
            sub=bind_rows(sub, readinfo)
        }
    }
    print(chrblock$chrname[1])
    return(sub)
}
cluster_copy(cluster, 'checkfilt')

read_filter <- function(i, cluster) { 
    ##read in motifs info
    motiffile=file.path(srcdir, paste0('barcodes', i, '.txt'))
    motifinfo=read_tsv(motiffile, col_names=FALSE)
    bc_cols=c('readname', 'chrname', motifinfo$X1)

    ##read in data
    barcodefile=file.path(datadir, paste0('20190809_zymo_control_barcodes', i,'.txt'))
    bcinfo=read_tsv(barcodefile, col_names=bc_cols, na=c('NA', 'None'))    
    countsfile=file.path(datadir, paste0('20190809_zymo_control_motifcounts', i,'.txt'))
    bccounts=read_tsv(countsfile, col_names=bc_cols, na=c('NA', 'None'))

    ##count NA and filter motif selections
    nacount=colSums(is.na(bcinfo)/dim(bcinfo)[1])
    lowna=nacount[nacount<.2]
    keepmotifs=names(lowna)
    bccounts=bccounts %>%
        select(all_of(keepmotifs))
    bcinfo=bcinfo %>%
        select(all_of(keepmotifs))

    ##filter out non-bacterial genome reads and any with less than 3 motifs
    countsfilter=bccounts %>%
        filter(!grepl('tig', chrname, fixed=TRUE)) %>%
        filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
        filter(across(c(-readname, -chrname), ~ .x>=3))
    cluster_copy(cluster, 'countsfilter')
    
    bcfilt=bcinfo %>%
        filter(!grepl('tig', chrname, fixed=TRUE)) %>%
        filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
        filter(complete.cases(.)) %>%
        group_by(chrname) %>%
        partition(cluster)
    bcfilt=bcfilt %>%
        do(checkfilt(., 10000, countsfilter)) %>%
        collect() %>%
        ungroup()
    
    return(bcfilt)
}

pcastuff <- function(bcfilt) {
    ##filter out cols
    bcdata=bcfilt %>%
        select(-readname, -chrname)
    
    ##scale and pca
    bcscaled=as_tibble(apply(bcdata, 2, scale))
    bcpca=prcomp(bcscaled)
    return(bcpca)
}


plotstuff <- function(bcfilt, pcasub, name, xlim=NULL, ylim=NULL) {
    ##umap of bcfilt
    bcdata=bcfilt %>%
        select(-readname, -chrname)
    bcumap=umap(bcdata)

    umapfile=file.path(dbxdir, paste0(name, '_umap.pdf'))
    pdf(umapfile, h=7, w=13)
    plot=plot_umap(bcumap$layout, bcfilt$chrname)
    print(plot)
    sep=plot + facet_wrap(~label)
    print(sep)
    dev.off()

    ##seurat workflow
    pcaumap=umap(pcasub)
    seuratfile=file.path(dbxdir, paste0(name, '_seurat.pdf'))
    pdf(seuratfile, h=7, w=13)
    plot=plot_umap(pcaumap$layout, bcfilt$chrname)
    print(plot)
    sep=plot + facet_wrap(~label)
    print(sep)
    if (!missing(xlim)) {
        plot=plot_umap(pcaumap$layout, bcfilt$chrname, xlim, ylim)
        print(plot)
        sep=plot + facet_wrap(~label)
        print(sep)
    }
    dev.off()
}



##top 10 most common motifs as barcodes
bcfilt=read_filter(conds[1], cluster)
bcfilt=bcfilt[complete.cases(bcfilt),]
bcpca=pcastuff(bcfilt)
pcasub=bcpca$x[,1:5]
plotstuff(bcfilt, pcasub, 'common10', c(-25, 25), c(-25,25))

bcfilt2=read_filter(conds[2], cluster)
bcpca2=pcastuff(bcfilt2)
pcasub2=bcpca2$x[,1:5]
plotstuff(bcfilt2, pcasub2, 'common15', c(-25, 25), c(-25,25))

bcfilt3=read_filter(conds[3], cluster)
bcpca3=pcastuff(bcfilt3)
pcasub3=bcpca3$x[,1:5]
plotstuff(bcfilt3, pcasub3, 'common20', c(-25, 25), c(-25,25))


