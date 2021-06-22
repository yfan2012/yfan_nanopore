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


read_filter <- function(i, cluster, minoccur) { 
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
        filter(across(c(-readname, -chrname), ~ .x>=minoccur))
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

read_plasmid <- function(i) {
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
        filter(grepl('plasmid', chrname, fixed=TRUE)) %>%
        filter(across(c(-readname, -chrname), ~ .x>=2))
    
    plas=bcinfo %>%
        filter(grepl('plasmid', chrname, fixed=TRUE)) %>%
        filter(complete.cases(.)) %>%
        filter(chrname %in% countsfilter$chrname)

    return(plas)
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
bcfilt=read_filter(conds[1], cluster, 3)
bcfilt=bcfilt[complete.cases(bcfilt),]
bcpca=pcastuff(bcfilt)
pcasub=bcpca$x[,1:5]
plotstuff(bcfilt, pcasub, 'common10', c(-25, 25), c(-25,25))

bcfilt2=read_filter(conds[2], cluster, 3)
bcpca2=pcastuff(bcfilt2)
pcasub2=bcpca2$x[,1:5]
plotstuff(bcfilt2, pcasub2, 'common15', c(-25, 25), c(-25,25))

bcfilt3=read_filter(conds[3], cluster, 3)
bcpca3=pcastuff(bcfilt3)
pcasub3=bcpca3$x[,1:5]
plotstuff(bcfilt3, pcasub3, 'common20', c(-25, 25), c(-25,25))




##plasmid analysis - take min occurence requirement to 2 or else you don't get any plasmid representation
bcfilt=read_filter(conds[2], cluster, 2)
bcfilt=bcfilt[complete.cases(bcfilt),]
plasfilt=read_plasmid(conds[2])
allfilt=bind_rows(bcfilt, plasfilt)
allpca=pcastuff(allfilt)
pcasub=allpca$x[,1:5]
plotstuff(allfilt, pcasub, 'plasmid')

alldata=allfilt %>%
    select(-chrname, -readname)
allumap=umap(alldata)
embed=tibble(x=allumap$layout[,1],
             y=allumap$layout[,2],
             label=allfilt$chrname) %>%
    mutate(colour=case_when(!grepl('plasmid', label, fixed=TRUE) ~ label, TRUE ~ 'plasmid')) %>%
    mutate(shape=case_when(grepl('plasmid', label, fixed=TRUE) ~ label, TRUE ~ 'chr'))
embedplas=embed %>%
    filter(colour=='plasmid') %>%
    select(-label)
embedchr=embed %>%
    filter(shape=='chr')


mycolors=c(brewer.pal(8, 'Set2'), '#000000')
myshapes=c(3,18)
plasfile=file.path(dbxdir, 'plasmid_shown.pdf')
pdf(plasfile, h=9, w=13)
plot=ggplot(embedchr, aes(x=x, y=y, colour=colour)) +
    geom_point(alpha=.2, size=.1) +
    scale_colour_manual(values=mycolors) +
    theme_bw()
mainplot=plot +
    geom_point(data=embedplas, aes(x=x, y=y, shape=shape), inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes) +
    theme(legend.position = 'none')
print(mainplot)
sep=plot +
    facet_wrap(~label) +
    geom_point(data=embedplas, aes(x=x, y=y, shape=shape), size=.5, alpha=.4, inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes) +
    theme(legend.position = "none")
print(sep)
dev.off()

