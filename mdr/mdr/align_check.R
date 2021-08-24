library(tidyverse)

##look at how useful the new reference is

datadir='/mithril/Data/Nanopore/projects/methbin/mdr'
prefix='200708_mdr_stool16native'
paffile=file.path(datadir, 'align', paste0(prefix, '.paf'))

pafcols=c('qname', 'qlen', 'qstart', 'qend', 'strand', 'rname', 'rlen', 'rstart', 'rend', 'num_match', 'alen', 'mapq', 'tp', 'cm', 's1', 's2', 'dv', 'rl')
paf=read_tsv(paffile, col_names=pafcols)

paffilt=paf %>%
    filter(mapq>45) %>%
    filter(alen>3000) %>%
    filter(num_match/alen>.9) %>%
    select(-tp, -cm, -s1, -s2, -dv, -rl)


