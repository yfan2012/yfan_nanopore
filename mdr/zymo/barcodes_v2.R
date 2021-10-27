library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')
projdir='/mithril/Data/Nanopore/projects/methbin/zymo'
datadir=file.path(projdir, 'barcode_v2', '20190809_zymo_control')
srcdir='~/Code/yfan_nanopore/mdr/rebase'
dbxdir='~/gdrive/mdr/zymo'

i=20
cluster_copy(cluster, 'checkfilt')


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
    filter(across(c(-readname, -chrname), ~ .x>=2))
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

plasfilt=bcinfo %>%
    filter(grepl('plasmid', chrname, fixed=TRUE)) %>%
    filter(complete.cases(.))

allfilt=bind_rows(bcfilt, plasfilt)

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
plasfile=file.path(dbxdir, 'plasmid_shown_v2.pdf')
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


####plasmid read classification by euclidean distance
cluster_copy(cluster, 'get_read_distances')
cluster_copy(cluster, 'classify_by_top_reads')
cluster_copy(cluster, 'bcfilt')

plasvote=plasfilt %>%
    rowwise() %>%
    partition(cluster) %>%
    do(classify_by_top_reads(., bcfilt, 50)) %>%
    collect() %>%
    ungroup()

voteinfocsv=file.path(projdir, 'read_classification', 'voteinfo_distance_refbased_v2.csv')
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

plascountspdf=file.path(dbxdir, 'classify_plasmid_nearest_dist_refbased_v2.pdf')
pdf(plascountspdf, w=15, h=7)
plot=ggplot(votecounts, aes(x=shortclass, y=frac, colour=shortclass, fill=shortclass, alpha=.5)) +
    geom_bar(stat='identity') +
    facet_wrap(~chrname) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    ylim(0,1) +
    theme_bw()
print(plot)
dev.off()






