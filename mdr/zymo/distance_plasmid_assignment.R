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


tigsbulk=scaledinfo %>%
    filter(type=='chr') %>%
    select(-type) %>%
    group_by(chrname) %>%
    summarise_if(is.numeric, mean)

source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')
plasdists=scaledplas %>%
    rowwise() %>%
    do(classify_plasmid_reads_distance(., tigsbulk))


##there is a more elegant way to do this
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


plascountspdf=file.path(dbxdir, 'classify_plasmid_majority_dist.pdf')
pdf(plascountspdf, w=15, h=7)
plot=ggplot(plascounts, aes(x=name, colour=name, fill=name, alpha=.5)) +
    geom_bar() +
    facet_wrap(~chrname) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()
