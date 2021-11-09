library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')

cluster=new_cluster(7)
cluster_library(cluster, 'tidyverse')
projdir='/mithril/Data/Nanopore/projects/methbin/zymo'
datadir=file.path(projdir, 'barcode')
srcdir='~/Code/yfan_nanopore/mdr/rebase'
dbxdir='~/gdrive/mdr/zymo'

conds=c('10', '15', '20')
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



####most common motifs as barcodes
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




####plasmid analysis
##take min occurence requirement to 2 or else you don't get any plasmid representation
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
pdf(plasfile, h=9, w=19)
plot=ggplot(embedchr, aes(x=x, y=y, colour=colour)) +
    geom_point(alpha=.2, size=.1) +
    scale_colour_manual(values=mycolors) +
    theme_bw()
mainplot=plot +
    geom_point(data=embedplas, aes(x=x, y=y, shape=shape), inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes)
print(mainplot)
sep=plot +
    facet_wrap(~label) +
    geom_point(data=embedplas, aes(x=x, y=y, shape=shape), size=.5, alpha=.4, inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes) +
    theme(legend.position = "none")
print(sep)
plot=ggplot(embedchr, aes(x=x, y=y, colour=colour)) +
    geom_point(alpha=.2, size=.1) +
    scale_colour_manual(values=mycolors) +
    theme_bw()
print(plot)
dev.off()

##contig analysis
tiginfo=allfilt %>%
    select(-readname) %>%
    group_by(chrname) %>%
    summarise_each(funs(mean))
clusttig=as.matrix(tiginfo %>% select(-chrname))
rownames(clusttig)=tiginfo$chrname
tigscaled=scale(as.data.frame(tiginfo[,2:12]))
rownames(tigscaled)=tiginfo$chrname

heatmappdf=file.path(dbxdir, 'heatmap_ref.pdf')
pdf(heatmappdf, h=7, w=14)
hm=heatmap(tigscaled, scale='none')
noscale=heatmap(clusttig, scale='none')
distplot=plot(hclust(dist(tigscaled)))
tiginfodist=plot(hclust(dist(clusttig)))
print(hm)
print(noscale)
print(distplot)
print(tiginfodist)
dev.off()




####plasmid read classification by umap neighbors
embedchr=embedchr %>%
    select(-shape, -colour)
embedplas=embedplas %>%
    select(-colour) %>%
    rename(label = shape)


classinfo=classify_umap_neighbors(embedplas, embedchr, 50)
ecoliclass=classinfo %>%
    filter(label=='Escherichia_coli_plasmid') %>%
    gather(key='bacteria', value='count', -x, -y, -label)
staphclass=classinfo %>%
    filter(label=='Staphylococcus_aureus_plasmid1') %>%
    gather(key='bacteria', value='count', -x, -y, -label)

votekey=names(classinfo[4:10])
classinfo$vote=votekey[max.col(classinfo[,4:10])]
majority=classinfo %>%
    select(label, vote) %>%
    rowwise() %>%
    mutate(short=strsplit(vote, '_', fixed=TRUE)[[1]][1]) %>%
    mutate(shorter=strsplit(short, '.', fixed=TRUE)[[1]][1])

majoritycount=majority %>%
    group_by(shorter, label) %>%
    summarise(count=n()) %>%
    ungroup() %>%
    group_by(label) %>%
    mutate(frac=count/sum(count))


classifypdf=file.path(dbxdir, 'classify_plasmid_umap_refbased.pdf')
pdf(classifypdf, w=15, h=7)
plot=ggplot(majoritycount, aes(x=shorter, y=frac, colour=shorter, fill=shorter, alpha=.5)) +
    geom_bar(stat='identity') +
    facet_wrap(~label) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()


####plasmid read classification by euclidean distance
cluster=new_cluster(36)
cluster_library(cluster, 'tidyverse')
cluster_copy(cluster, 'get_read_distances')
cluster_copy(cluster, 'classify_by_top_reads')
cluster_copy(cluster, 'bcfilt')

plasvote=plasfilt %>%
    rowwise() %>%
    partition(cluster) %>%
    do(classify_by_top_reads(., bcfilt, 50)) %>%
    collect() %>%
    ungroup()

voteinfocsv=file.path(projdir, 'read_classification', 'voteinfo_distance_refbased.csv')
write_csv(plasvote, voteinfocsv)
plasvote=read_csv(voteinfocsv)


voteclass=tibble(chrname=as.character(),
                 class=as.character())

for (i in 1:dim(plasvote)[1]) {
    readdata=plasvote[i,]
    nearest=colnames(readdata)[1:7][which(readdata[1:7]==max(readdata[1:7]))]

    readclass=tibble(chrname=readdata$chrname, class=nearest)
    voteclass=bind_rows(voteclass, readclass)
}


votecounts=voteclass %>%
    group_by(chrname, class) %>%
    summarise(counts=n()) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(abrevclass=strsplit(class, split='_', fixed=TRUE)[[1]][1]) %>%
    mutate(shortclass=strsplit(abrevclass, split='.', fixed=TRUE)[[1]][1]) %>%
    select(-abrevclass, -class) %>%
    group_by(chrname) %>%
    mutate(frac=counts/sum(counts))

plascountspdf=file.path(dbxdir, 'classify_plasmid_nearest_dist_refbased.pdf')
pdf(plascountspdf, w=15, h=7)
plot=ggplot(votecounts, aes(x=shortclass, y=frac, colour=shortclass, fill=shortclass, alpha=.5)) +
    geom_bar(stat='identity') +
    facet_wrap(~chrname) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()


refbulk=bcfilt %>%
    group_by(chrname) %>%
    summarise_if(is.numeric, mean)

plasdists=plasfilt %>%
    rowwise() %>%
    do(classify_plasmid_reads_distance(.,refbulk))

plasclass=tibble(chrname=as.character(),
                 plasclass=as.character())
for (i in 1:dim(plasdists)[1]) {
    readdata=plasdists[i,]
    nearest=colnames(readdata)[-1][which(readdata[-1]==min(readdata[-1]))]

    readclass=tibble(chrname=readdata$chrname, plasclass=nearest)
    plasclass=bind_rows(plasclass, tibble(chrname=readdata$chrname, plasclass=nearest))
}

plascounts=plasclass %>%
    rowwise() %>%
    mutate(tally=strsplit(chrname, '_', fixed=TRUE)[[1]][1]==strsplit(plasclass, '_', fixed=TRUE)[[1]][1]) %>%
    mutate(short=strsplit(plasclass, '_', fixed=TRUE)[[1]][1]) %>%
    mutate(name=strsplit(short, '.', fixed=TRUE)[[1]][1])

plasvotes=plascounts %>%
    group_by(chrname, name) %>%
    summarise(counts=n()) %>%
    ungroup() %>%
    group_by(chrname) %>%
    mutate(frac=counts/sum(counts))

plascountspdf=file.path(dbxdir, 'classify_plasmid_pseudobulk_dist_refbased.pdf')
pdf(plascountspdf, w=15, h=7)
plot=ggplot(plasvotes, aes(x=name, y=frac, colour=name, fill=name, alpha=.5)) +
    geom_bar(stat='identity') +
    facet_wrap(~chrname) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()
