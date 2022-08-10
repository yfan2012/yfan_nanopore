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
library(cluster)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')


projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'

chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

invec=c('CAGAG','CCWGG', 'CMTCGAKG','CTCCAG', 'CTKVAG', 'GATC', 'GCGC', 'GCWGC', 'GGCC', 'GGNNCC', 'GGWCC', 'TCCGGA', 'GCCGGC', 'RGCGCY')
methfreq=readRDS(file.path(datadir, 'methfreq.rds')) %>%
    rowwise() %>%
    filter(motif %in% invec)
    

fullfreqs=methfreq %>%
    spread(key=motif, value=freq)
fullnacount=colSums(is.na(fullfreqs))/dim(fullfreqs)[1]
fullcounts=tibble(motifs=names(fullnacount), frac=1-fullnacount)
fullcountcsv=file.path(dbxdir, 'motifcounts.csv')
write_csv(fullcounts, fullcountcsv)


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
    ##diana %>%
    hclust %>%
    as.dendrogram
dendfile=file.path(dbxdir, 'plaindend.hclust.tree')
write.tree(dend %>% as.phylo, dendfile)

##try div instead of agg
ddend=matchrominfo %>%
    scale %>%
    dist %>%
    diana %>%
    as.dendrogram
ddendfile=file.path(dbxdir, 'plaindend.diana.tree')
write.tree(ddend %>% as.phylo, ddendfile)


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
    mutate(type=case_when(label %in% unique(plasnear$chroms) ~ label,
                                                       TRUE ~ ''))
dendlabeled=dend %>%
    set('leaves_pch', 15) %>%
    set('leaves_cex', 2) %>%
    set('leaves_col', labelinfo$color) %>%
    set('labels', labelinfo$type)

labtreepdf=file.path(dbxdir, 'tree_labeled.pdf')
pdf(labtreepdf, h=35, w=8)
par(mar = c(3, 4, 2, 15) + 0.1)
plot(dendlabeled, horiz=TRUE)
dev.off()

heatmappdf=file.path(dbxdir, 'heatmap.pdf')
pdf(heatmappdf, h=18, w=18)
heatmap(matchrominfo)
dev.off()



####incorporating kraken analysis
krakentigfile=file.path(datadir, 'tigbin_species_inclusive_kraken.tsv')
krakentig=read_tsv(krakentigfile)

krakenlabel=labelinfo %>%
    rowwise() %>%
    mutate(tigspecies=krakentig$tigkraken[krakentig$tig==label]) %>%
    mutate(specieslab=strsplit(tigspecies, ' (', fixed=TRUE)[[1]][1]) %>%
    mutate(fulllab=paste0(c(label, specieslab), collapse=' '))

dendspecies=dend %>%
    set('leaves_pch', 15) %>%
    set('leaves_cex', 2) %>%
    set('leaves_col', krakenlabel$color) %>%
    set('labels', krakenlabel$fulllab)

speciestreepdf=file.path(dbxdir, 'tree_labeled_species.pdf')
pdf(speciestreepdf, h=35, w=12)
par(mar = c(3, 4, 2, 23) + 0.1)
plot(dendspecies, horiz=TRUE)
dev.off()

