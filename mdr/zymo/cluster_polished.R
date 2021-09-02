library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/mdr/cluster_functions.R')


datadir='/mithril/Data/Nanopore/projects/methbin/zymo'
dbxdir='~/gdrive/mdr/mdr'
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
    filter(across(c(-readname, -chrname), ~.x>=3))

bcfilt=bcinfo %>%
    filter(readname %in% bccounts$readname) %>%
    filter(chrname %in% key$tig) %>%
    rowwise() %>%
    mutate(chrname=key$species[key$tig==chrname])
 
plasinfo=bcfilt %>%
    filter(chrname %in% plasnames)
