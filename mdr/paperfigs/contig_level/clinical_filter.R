library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

##use coverage to filter out list of potentially misassembled contigs
##don't want to lose too many plasmid contigs, so take that into consideration 

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_asm'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'








####try using mean deviation from average
covdev=cov %>%
    group_by(tig) %>%
    mutate(avgcov=mean(cov)) %>%
    mutate(dev=abs(cov-avgcov)/avgcov) %>%
    summarise(avgdev=mean(dev)) %>%
    arrange(-avgdev) %>%
    mutate(rank=dense_rank(desc(avgdev))) %>%
    rowwise() %>%
    mutate(len=covpeaks$len[covpeaks$tig==tig])




####just try a q-q plot??
##qq plot and use some kind of residual 
covqq=covfreq %>%
    do(qqmse(.)) %>%
    filter(tig %in% label_order) %>%
    arrange(-mse)
##doesn't look as compeling as my weird thing looking for bimodality 
##tails are messing with me?
