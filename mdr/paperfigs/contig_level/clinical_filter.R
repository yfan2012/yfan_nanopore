library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

##use coverage to filter out list of potentially misassembled contigs
##don't want to lose too many plasmid contigs, so take that into consideration 

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_asm'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'






####coverage hist analysis
covfile=file.path(projdir, 'mdr', 'align', paste0(prefix, 'polished.sorted.cov'))
cov_cols=c('tig', 'pos', 'cov')
cov=read_tsv(covfile, col_names=cov_cols) %>%
    group_by(tig) %>%
    mutate(covfrac=cov/max(cov))


##get coverage spectrum
covfreq=cov %>%
    group_by(tig, cov) %>%
    summarise(freq=n()) %>%
    mutate(normfreq=freq/max(freq)) %>%
    mutate(normcov=cov/max(cov))
    
covpeaks=covfreq %>%
    do(findpeaks(.)) %>%
    mutate(composite=sum(h*t*s)) %>%
    arrange(-composite) %>%
    ungroup() %>%
    mutate(rank=dense_rank(desc(composite)))



####get list of plasmid tigs
tigplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.plasmidfinder.tsv')
plas_cols=c('file', 'seq', 'start', 'end', 'strand', 'gene', 'coverage', 'covmap', 'gaps', 'covfrac', 'ident', 'db',
            'acc', 'prod', 'res')
tigplas=read_tsv(tigplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc) %>%
    rowwise() %>%
    mutate(rank=covpeaks$rank[covpeaks$tig==seq]) %>%
    arrange(rank)








####just try a q-q plot??
##qq plot and use some kind of residual 
covqq=covfreq %>%
    do(qqmse(.)) %>%
    filter(tig %in% label_order) %>%
    arrange(-mse)
##doesn't look as compeling as my weird thing looking for bimodality 
##tails are messing with me?
