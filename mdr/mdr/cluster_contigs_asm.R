library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/mdr/cluster_functions.R')

##cluster=new_cluster(6)
##cluster_library(cluster, 'tidyverse')
##cluster_copy(cluster, 'checkfilt')

datadir='/mithril/Data/Nanopore/projects/methbin/mdr'
dbxdir='~/gdrive/mdr/mdr'
prefix='200708_mdr_stool16native'
prefixasm=paste0(prefix, '_asm')

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)

tigidfile=file.path(datadir, 'contig_id', 'blast_tigid.csv')
tigid=read_csv(tigidfile)
tigdups=tigid %>%
    group_by(tig) %>%
    mutate(counts=n()) %>%
    filter(counts>1) %>%
    select(c(-refmax, -counts))
tigiduniq=tigid %>%
    group_by(tig) %>%
    mutate(counts=n()) %>%
    filter(counts==1) %>%
    select(c(-refmax, -counts))

barcodedatafile=file.path(datadir, 'barcode', prefixasm, paste0(prefix, '_barcodes50.txt'))
fullbcinfo=read_tsv(barcodedatafile, col_names=bc_cols, na=c('None'))
barcodecountsfile=file.path(datadir, 'barcode', prefixasm, paste0(prefix, '_barcodes50_motifcounts.txt'))
fullbccounts=read_tsv(barcodecountsfile, col_names=bc_cols)

nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
lowna=names(nacount[nacount<.2])

bcinfo=fullbcinfo %>%
    select(all_of(lowna))
bccounts=fullbccounts %>%
    select(all_of(lowna)) %>%
    filter(complete.cases(.)) %>%
    filter(across(c(-readname, -chrname), ~.x>=5))

bcfilt=bcinfo %>%
    filter(complete.cases(.)) %>%
    filter(readname %in% bccounts$readname)

tiginfo=bcfilt %>%
    group_by(chrname) %>%
    summarise_if(is.numeric, mean) %>%
    mutate(genus=case_when(chrname %in% tigiduniq$tig ~ tigiduniq$gmax[which(chrname==tigiduniq$tig)],
                           TRUE ~ 'unknown'))
tigdata=tiginfo %>%
    select(-chrname)

tigumap=umap(tigdata)
tigembed=tibble(x=tigumap$layout[,1],
                y=tigumap$layout[,2],
                chrname=tiginfo$chrname)



tigplot=file.path(dbxdir, 'contig_umap.pdf')
pdf(tigplot, h=9, w=13)
##umap of contigs, no colours
plot=ggplot(tigembed, aes(x=x, y=y)) +
    geom_point(alpha=.4, size=.5) +
    theme_bw()
print(plot)
dev.off()
