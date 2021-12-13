library(tidyverse)
library(multidplyr)
library(Biostrings)

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

dbxdir='~/gdrive/mdr/etc'
projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'etc/CG_model')

mbedcols=c('chr', 'start', 'end', 'name', 'score', 'strand', 'cstart', 'cend','color', 'cov', 'meth')
methcallsfile=file.path(datadir, 'modified_bases.5mC.bed')


####read motif info
cmethfile=file.path(projdir, 'zymo/truth/bisulfite/zymo_cmeth.csv')
cmeth=read_csv(cmethfile)
amethfile=file.path(projdir, 'zymo/truth/pacbio/zymo_ameth.csv')
ameth=read_csv(amethfile)
methinfo=bind_rows(cmeth,ameth) %>%
    group_by(motif, pos) %>%
    summarise(num=n())

nametochrfile='~/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
chrlabels=read_table2(nametochrfile, col_names=c('chr', 'label'))

meth=read_tsv(methcallsfile, col_names=mbedcols) %>%
    filter(!grepl('tig', chr, fixed=TRUE))



####bisulf stuff
keyfile='~/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
keycols=c('chrom', 'label')
key=read_table2(keyfile, col_names=keycols)

labels=c('bsubtilis', 'ecoli', 'efaecalis', 'lmonocytogenes', 'paeruginosa', 'saureus', 'senterica')
cxcols=c('chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'seq')
allcx=NULL
for (label in labels) {
    cxfile=file.path(projdir, 'zymo/truth/bisulfite/bismark', label, paste0(label, '_1_bismark_bt2_pe.CX_report.txt'))
    chrlist=key$chrom[key$label==label]
    cx=read_tsv(cxfile, col_names=cxcols) %>%
        rowwise() %>%
        filter(chr %in% chrlist) %>%
        mutate(total=meth+unmeth)
    allcx=bind_rows(allcx, cx)
}

methbisulf=allcx %>%
    mutate(bisulf=meth/total) 




####filter out appropriate 
reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
ref=readDNAStringSet(reffile, format='fasta')
chrs=names(ref)[!grepl('tig', names(ref), fixed=TRUE)]
methcompare=NULL
for (i in chrs) {
    print(i)
    chromname=strsplit(i, split=' ', fixed=TRUE)[[1]][1]
    starts=start(vmatchPattern('GCCGGC', ref[i])[[i]])
    bisulfpos=c(starts+2, starts+3)
    bisulf=methbisulf %>%
        filter(chr==chromname) %>%
        rowwise() %>%
        filter(pos %in% bisulfpos) %>%
        mutate(motif='GCCGGC')
    mega=meth %>%
        filter(chr==chromname) %>%
        filter(end %in% bisulfpos) %>%
        mutate(pos=end)
    allmeth=full_join(bisulf, mega, by=c('pos')) %>%
        select(chr.x, pos, bisulf, meth.y) %>%
        rename(chr.x='chr') %>%
        rename(meth.y='mega')

    methcompare=bind_rows(methcompare, allmeth)
}

    
boxpdf=file.path(dbxdir, 'mega_boxes_bed.pdf')
pdf(boxpdf, h=9, w=13)
box=ggplot(methcompare, aes(x=chr, y=mega, colour=chr, fill=chr, alpha=.2)) +
    geom_boxplot() +
    ggtitle('GCCGGC') +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(box)
dev.off()


corpdf=file.path(dbxdir, 'bisulf_cor_bed.pdf')
pdf(corpdf, h=9, w=13)
box=ggplot(methcompare, aes(x=bisulf, y=mega, colour=chr, alpha=.2)) +
    geom_point() +
    ggtitle('GCCGGC') +
    theme_bw()
print(box)
dev.off()

