library(tidyverse)
source('~/Code/yfan_nanopore/mdr/disco/contig_id_functions.R')

prefix='20190809_zymo_control_polished'
datadir=file.path('/mithril/Data/Nanopore/projects/methbin/zymo/barcode', prefix)


blasttsv=file.path(datadir, paste0(prefix, '_blast.tsv'))
blastcols=c('tig', 'ref', 'ident', 'alen', 'mismatch', 'gaps', 'tigstart', 'tigend', 'refstart', 'refend', 'eval', 'bit')

blast=read_tsv(blasttsv, col_names=blastcols, comment='#') %>%
    select(-c(mismatch, gaps, refstart, refend, eval, bit)) %>%
    filter(alen>5000 & ident>98)

merged=blast %>%
    group_by(tig, ref) %>%
    do(merge_overlaps(.)) %>%
    mutate(org=case_when(!grepl('tig', ref, fixed=TRUE) ~ ref, TRUE ~ 'yeast'))    

key=merged %>%
    group_by(tig) %>%
    summarise(species=case_when(length(table(org))==1 ~ org[1], TRUE ~ 'mixed')) %>%
    filter(species!='yeast')


