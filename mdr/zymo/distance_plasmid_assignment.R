library(tidyverse)
library(umap)
library(lazyeval)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo'
dbxdir='~/gdrive/mdr/zymo'
prefix='20190809_zymo_control_polished'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)
keyfile=file.path(datadir, 'medaka', 'consensus_key.csv')
key=read_csv(keyfile)


barcodedatafile=file.path(datadir, 'barcode', prefix, '20190809_zymo_control_contigs_barcodes15.txt')
fullbcinfo=read_tsv(barcodedatafile, col_names=bc_cols, na=c('None'))
countsfile=file.path(datadir, 'barcode', prefix, '20190809_zymo_control_barcodes15_motifcounts.txt')
fullbccounts=read_tsv(countsfile, col_name=bc_cols, na=c('None'))

plasnames=c('Staphylococcus_aureus_plasmid1' ,'Escherichia_coli_plasmid')
plastigs=key$tig[key$species %in% plasnames]
plasbcinfo=fullbcinfo %>%
    filter(chrname %in% plastigs)

nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
lowna=nacount[nacount<.2]
plasnacount=colSums(is.na(plasbcinfo)/dim(plasbcinfo)[1])
plaslowna=plasnacount[plasnacount<.2]
keepmotifs=intersect(names(lowna), names(plaslowna))


bcinfo=fullbcinfo %>%
    select(all_of(keepmotifs))
bccounts=fullbccounts %>%
    select(all_of(keepmotifs)) %>%
    filter(complete.cases(.)) %>%
    filter(across(c(-readname, -chrname), ~.x>=5)) %>%
    filter(chrname %in% key$tig) %>%
    rowwise() %>%
    mutate(chrname=key$species[key$tig==chrname])

allbcfilt=bcinfo %>%
    filter(complete.cases(.)) %>%
    filter(readname %in% bccounts$readname) %>%
    filter(chrname %in% key$tig) %>%
    rowwise() %>%
    mutate(chrname=key$species[key$tig==chrname])

plasinfo=allbcfilt %>%
    filter(chrname %in% plasnames) %>%
    mutate(type='plas')

##standardize
bcfilt=allbcfilt %>%
    filter(!(chrname %in% plasnames)) %>%
    group_by(chrname) %>%
    do(checkfilt(., 10000, bccounts)) %>%
    ungroup() %>%
    mutate(type='chr')

allinfo=bind_rows(bcfilt, plasinfo)
scaledinfo=allinfo %>%
    mutate_if(is.numeric, ~(scale(.) %>% as.vector))
scaledplas=scaledinfo %>%
    filter(type=='plas') %>%
    select(-type)




##classify using min distance from pseudobulk tigs
tigsbulk=scaledinfo %>%
    filter(type=='chr') %>%
    select(-type) %>%
    group_by(chrname) %>%
    summarise_if(is.numeric, mean)

plasdists=scaledplas %>%
    rowwise() %>%
    do(classify_plasmid_reads_distance(., tigsbulk))


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
    

plascountspdf=file.path(dbxdir, 'classify_plasmid_pseudobulk_dist.pdf')
pdf(plascountspdf, w=15, h=7)
plot=ggplot(plasvotes, aes(x=name, y=frac, colour=name, fill=name, alpha=.5)) +
    geom_bar(stat='identity') +
    facet_wrap(~chrname) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()



##classify by using 50 closest reads by euclidean distance
cluster=new_cluster(36)
cluster_library(cluster, 'tidyverse')
cluster_copy(cluster, 'get_read_distances')
cluster_copy(cluster, 'classify_by_top_reads')
cluster_copy(cluster, 'bcfilt')

plasvote=plasinfo %>%
    rowwise() %>%
    partition(cluster)

voteinfo=plasvote %>%
    do(classify_by_top_reads(., bcfilt, 50))
voteresults=voteinfo %>%
    collect() %>%
    ungroup()

voteinfocsv=file.path(datadir, 'read_classification', 'voteinfo_distance_tigbased.csv')
write_csv(voteresults, voteinfocsv)


voteclass=tibble(chrname=as.character(),
                 class=as.character())

for (i in 1:dim(voteresults)[1]) {
    readdata=voteresults[i,]
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
    select(-abrevclass, -class)
votecounts=bind_rows(votecounts, tibble(chrname='Escherichia_coli_plasmid', counts=0, shortclass='Listeria'))
votecounts=votecounts %>%
    group_by(chrname) %>%
    mutate(frac=counts/sum(counts))

plascountspdf=file.path(dbxdir, 'classify_plasmid_nearest_dist.pdf')
pdf(plascountspdf, w=15, h=7)
plot=ggplot(votecounts, aes(x=shortclass, y=frac, colour=shortclass, fill=shortclass, alpha=.5)) +
    geom_bar(stat='identity') +
    facet_wrap(~chrname) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()
