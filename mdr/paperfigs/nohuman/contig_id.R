library(tidyverse)

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'


####assign contigs to bins based on the mummer
mumfile=file.path(datadir, 'clin_mummer', 'asm_hiC.1coords')
mumcols=c('rstart', 'rend', 'qstart', 'qend', 'ralen', 'qalen','ident', 'rlen', 'qlen', 'rcov', 'qcov', 'rname', 'qname')
mum=read_tsv(mumfile, col_names=mumcols) %>%
    rowwise() %>%
    mutate(bin=strsplit(qname, '.', fixed=TRUE)[[1]][1]) %>%
    group_by(rname, bin, rlen) %>%
    summarise(total_rcov=sum(rcov)) %>%
    filter(total_rcov>20)

mumkey=mum %>%
    group_by(rname) %>%
    filter(n()==1)

mummultiple=mum %>%
    group_by(rname) %>%
    filter(n()>1) %>%
    summarise(all_bins=paste0(bin, collapse=','), all_rcov=paste0(total_rcov, collapse=','))

mumkeyfile=file.path(datadir, 'tigs2bins.tsv')
write_tsv(mumkey, mumkeyfile)
multikeyfile=file.path(datadir, 'tigs2bins_multi.tsv')
write_tsv(mummultiple, multikeyfile)