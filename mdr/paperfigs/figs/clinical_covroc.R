library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(fossil)
library(mclust)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/figs'


####methylation distance
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))

plotallroc=read_csv(file.path(dbxdir, 'plotallroc.csv')) %>%
    mutate(randtype='leaf') %>%
    filter(samp<=5)
plotallrands=read_csv(file.path(dbxdir, 'plotallrands.csv')) %>%
    mutate(randtype='tree') %>%
    filter(samp<=5) %>%
    filter(label=='rando')
rocplot=bind_rows(plotallroc, plotallrands) %>%
    mutate(samp=as.character(samp))

covrocpdf=file.path(dbxdir, 'roc.pdf')
pdf(covrocpdf, h=8, w=11)
seqplot=ggplot(rocplot %>% filter(label=='real'), aes(x=seqtogether, y=seqpure, colour=samp)) +
    geom_step() +
    geom_point(rocplot %>% filter(label=='rando'), mapping=aes(x=seqtogether, y=seqpure, colour=randtype, alpha=.02)) +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on amount of sequence') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(seqplot)
tigplot=ggplot(rocplot %>% filter(label=='real'), aes(x=numtogether, y=numpure, colour=samp)) +
    geom_step() +
    geom_point(rocplot %>% filter(label=='rando'), mapping=aes(x=numtogether, y=numpure, colour=randtype, alpha=.02)) +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on number of contigs') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(tigplot)
dev.off()

