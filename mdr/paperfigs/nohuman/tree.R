library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(fossil)
library(mclust)
library(ggtree)
library(treeio)
library(ape)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')


projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'

chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

invec=c('CAGAG','CCWGG', 'CMTCGAKG','CTCCAG', 'CTKVAG', 'GATC', 'GCGC', 'GCWGC', 'GGCC', 'GGNNCC', 'GGWCC', 'TCCGGA', 'GCCGGA', 'RGCGCY')
methfreq=readRDS(file.path(datadir, 'methfreq.rds')) %>%
    rowwise() %>%
    filter(motif %in% invec)
    

fullfreqs=methfreq %>%
    spread(key=motif, value=freq)
fullnacount=colSums(is.na(fullfreqs))/dim(fullfreqs)[1]
fullcounts=tibble(motifs=names(fullnacount), frac=1-fullnacount)

freqs=methfreq %>%
    rowwise() %>%
    filter(motif %in% invec) %>%
    spread(key=motif, value=freq)
nacount=colSums(is.na(freqs))/dim(freqs)[1]


nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])

methchroms=methfreq %>%
    filter(chrom %in% keepchroms)
chrominfo=methchroms %>%
    spread(key=motif, value=freq) %>%
    filter(chrom %in% chrombins$rname)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom
dend=matchrominfo %>%
    scale %>%
    dist %>%
    hclust %>%
    as.dendrogram
dendfile=file.path(dbxdir, 'plaindend.tree')
write.tree(dend %>% as.phylo, dendfile)

mycolors=read_csv('~/Code/yfan_nanopore/mdr/paperfigs/figs/colors.csv')
mycolors=rbind(mycolors, c('bin_10', '#000000'))
mycolors=rbind(mycolors, c('bin_18', '#FFFFFF'))

labelinfo=tibble(label=labels(dend)) %>%
    filter(label %in% chrombins$rname) %>%
    rowwise() %>%
    mutate(bins=chrombins$bin[chrombins$rname==label]) %>%
    mutate(color=mycolors$color[mycolors$bin==bins])

plaindend=dend %>%
    set('leaves_pch', 15) %>%
    set('leaves_cex', 2) %>%
    set('leaves_col', labelinfo$color) %>%
    set('labels', rep('', length(labelinfo$color)))

treepdf=file.path(dbxdir, 'tree.pdf')
pdf(treepdf, h=35, w=8)
par(mar = c(3, 4, 2, 15) + 0.1)
plot(plaindend, horiz=TRUE)
dev.off()



plasnearbinscsv=file.path(dbxdir, 'plasmid_nearest.csv')
plasnear=read_csv(plasnearbinscsv)
labelinfo=labelinfo %>%
    mutate(type=case <- when(label %in% unique(plasnear$chroms) ~ label,
                                                       TRUE ~ ''))
