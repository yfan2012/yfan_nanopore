library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')


projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'


methfile=file.path(datadir, 'methcalls_nohuman_perf.csv')
methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
meth=read_csv(methfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))

cluster_copy(cluster, 'findMethFreq')

methgrouped=meth %>%
    filter(sum(methnum+umethnum)>5) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methfreq=methgrouped %>%
    do(findMethFreq(.))  %>%
    collect() %>%
    summarise(freq=mean(methfrac))


methfreqrds=file.path(datadir, 'methfreq.rds')
saveRDS(methfreq, methfreqrds)
