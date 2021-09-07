library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/mdr/cluster_functions.R')


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
    filter(chrname %in% plasnames)

bcfilt=allbcfilt %>%
    filter(!(chrname %in% plasnames)) %>%
    group_by(chrname) %>%
    do(checkfilt(., 10000, bccounts)) %>%
    ungroup()

bcdata=rbind(bcfilt, plasinfo)



bcumap=umap(bcdata %>% select(-chrname, -readname))

embedinfo=tibble(x=bcumap$layout[,1],
             y=bcumap$layout[,2],
             label=bcdata$chrname)
embedchr=embedinfo %>%
    filter(!(label %in% plasnames))
embedplas=embedinfo %>%
    filter(label %in% plasnames)



mycolors=c(brewer.pal(8, 'Set2'), '#000000')
myshapes=c(3,18)
plasfile=file.path(dbxdir, 'plasmid_shown_polished.pdf')
pdf(plasfile, h=9, w=13)
plot=ggplot(embedchr, aes(x=x, y=y, colour=label)) +
    geom_point(alpha=.2, size=.1) +
    scale_colour_manual(values=mycolors) +
    theme_bw()
mainplot=plot +
    geom_point(data=embedplas, aes(x=x, y=y, shape=label), inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes) +
    theme(legend.position='none')
print(mainplot)
dev.off()

