library(tidyverse)

projdir='/mithril/Data/Nanopore/projects/methbin/zymo'
pbdir=file.path(projdir, 'truth/pacbio')
smrtdir=file.path(pbdir, 'smrtanalysis/016')

samps=c('016438','016439','016441','016442','016443','016444','016445','016446')

covcols=c('chrom', 'start', 'end', 'description', 'meancov', 'strand')


####assign each smrtanalysis number to a sample label
key=NULL
for (i in samps) {
    covfile=file.path(smrtdir, i, 'data/coverage.bed')
    cov=read_tsv(file=covfile, col_names=covcols, skip=1) %>%
        group_by(chrom) %>%
        summarise(total=mean(meancov)) %>%
        filter(total>20)
    cov$samp=i
    key=bind_rows(key, cov)
}

chrlistfile='~/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
chrlist=read_table2(chrlistfile, col_names=c('chr', 'label'))

labeledkey=key %>%
    rowwise() %>%
    mutate(chrom=strsplit(chrom, split=' ', fixed=TRUE)[[1]][1]) %>%
    mutate(label=chrlist$label[chrlist$chr==chrom]) %>%
    select(samp, label) %>%
    distinct()
labeledkey$label[8]='ecoli_ii'

pblabfile='~/Code/yfan_nanopore/mdr/zymo/truth/pblist_withlabels.txt'
write_tsv(labeledkey, pblabfile, col_names=FALSE)
