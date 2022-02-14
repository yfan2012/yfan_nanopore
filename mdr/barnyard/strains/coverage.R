library(tidyverse)

projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'barnyard/strains')

prefix='220131_mdr_barnyard_'
dbxdir='~/gdrive/mdr/paperfigs/contig_level'

covcols=c('chrom', 'pos', 'cov')
scovfile=file.path(datadir, 'align', paste0(prefix, 'st3294.genomecov'))
ecovfile=file.path(datadir, 'align', paste0(prefix, 'st3689.genomecov'))
scov=read_tsv(scovfile, col_names=covcols) %>%
    group_by(chrom) %>%
    summarise(avgcov=mean(cov))
ecov=read_tsv(ecovfile, col_names=covcols) %>%
    group_by(chrom) %>%
    summarise(avgcov=mean(cov))
