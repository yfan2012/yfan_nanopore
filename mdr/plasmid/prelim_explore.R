library(tidyverse)

datadir='/mithril/Data/Nanopore/projects/methbin/plasmid/blast'

####grab blast data and filter
filelist=list.files(datadir)
blastcols=c('query', 'subject', 'ident', 'alen', 'mm', 'gap', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bit')

fullblast=tibble()
for (i in filelist) {
    blastfile=file.path(datadir, i)
    blast=read_tsv(blastfile, col_names=blastcols, comment='#')
    fullblast=bind_rows(fullblast,blast)
}

rebasedir='/uru/Data/Nanopore/projects/mdr'
rebasecsv=file.path(rebasedir, 'refs', 'rebase_report.csv')
rebase=read_csv(rebasecsv) %>%
    filter(Gene=='M') %>%
    filter(Type=='II') %>%
    mutate(len=abs(Start-End))

filtblast=fullblast %>%
    rowwise() %>%
    mutate(mtaselen=min(rebase$len[which(rebase$SeqName==subject)])) %>%
    filter(ident>95) %>% 
    filter(alen/mtaselen>.95)


pyoblast=filtblast %>%
    select(query, subject, qstart, qend)

redblast=filtblast %>%
    group_by(query, subject) %>%
    summarise(alen=max(alen)) %>%
    ungroup() %>%
    mutate(motif=rebase$Specificity[rebase$SeqName==subject][1])


resfile=file.path(datadir, 'mtase_blast_filt.csv')
write_csv(redblast, resfile)

