library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'


####normal distance calcs
##methfile=file.path(datadir, 'clin_barocdes_methcalls.perf2.csv')
methfile=file.path(datadir, 'clin_barocdes_methcalls.perf3.csv')
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

##methfreqrds=file.path(datadir, 'clin_methfreq.rds')
methfreqrds=file.path(datadir, 'clin_methfreq3.rds')
saveRDS(methfreq, methfreqrds)



####coverage based methfrac
covcols=c('chrom', 'pos', 'cov')
covfile=file.path(projdir, 'mdr/align/200708_mdr_stool16native_asmpolished.perf.sorted.cov')
cov=read_tsv(covfile, covcols)

methcov=inner_join(meth, cov, by=c('chrom', 'pos'))

covgrouped=methcov %>%
    filter(sum(methnum+umethnum)>5) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
covloci=covgrouped %>%
    do(findMethFreq(.))  %>%
    collect()
covfreq=covloci %>%
    ungroup() %>%
    mutate(methfrac=methnum/cov) %>%
    group_by(motif, chrom) %>%
    summarise(freq=mean(methfrac), numloci=n())

covfreqrds=file.path(datadir, 'clin_methfreq3_cov.rds')
saveRDS(covfreq, covfreqrds)
    
