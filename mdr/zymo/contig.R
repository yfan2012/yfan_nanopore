library(tidyverse)

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/contig'
dbxdir='~/gdrive/mdr/zymo'
prefix='20190809_zymo_control'

methfile=file.path(datadir, paste0(prefix, '_contig15.txt'))
covfile=file.path(datadir, paste0(prefix, '_contig15_cov.txt'))

meth=read_tsv(methfile) %>%
    group_by(chrname) %>%
    summarise(across(everything(), sum)) %>%
    ungroup() %>%
    filter(!grepl('tig', chrname, fixed=TRUE))

cov=read_tsv(covfile) %>%
    group_by(chrname) %>%
    summarise(across(everything(), sum)) %>%
    ungroup() %>%
    filter(!grepl('tig', chrname, fixed=TRUE))

ratio=scale(as.data.frame(meth[,-1]/cov[,-1]))
rownames(ratio)=cov$chrname
ratio=ratio[1:9,] ##staph plasmid2 and plasmid3 cov is low

heatmappdf=file.path(dbxdir, 'heatmap_contig_calls.pdf')
pdf(heatmappdf, h=7, w=13)
hm=heatmap(ratio, scale='none')
print(hm)
dev.off()


