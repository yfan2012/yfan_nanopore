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

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/figs'


##chrom to bin info
chrombinsfile='/mithril/Data/Nanopore/projects/methbin/paperfigs/contig_level/tigs2bins.tsv'
chrombins=read_tsv(chrombinsfile)


####get tree from methylation distance
exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))
freqs=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec) %>%
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
plaindend=matchrominfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
dendfile=file.path(dbxdir, 'plaindend.tree')
write.tree(plaindend %>% as.phylo, dendfile)


labelinfo=tibble(label=labels(plaindend)) %>%
    filter(label %in% chrombins$rname) %>%
    rowwise() %>%
    mutate(bins=chrombins$bin[chrombins$rname==label])
numcolors=length(names(table(labelinfo$bins)))
colors=c(colorRampPalette(brewer.pal(8, 'Set2'))(numcolors+1))
colors=sample(colors[1:numcolors])
mycolors=tibble(color=colors,
                bin=names(table(labelinfo$bins)))

labelinfo=labelinfo %>%
    mutate(color=mycolors$color[mycolors$bin==bins])

nodePar=list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col='blue')
plaindend=plaindend %>%
    set('leaves_pch', 19) %>%
    set('leaves_cex', 2 ) %>%
    set('leaves_col', labelinfo$color) %>%
    set('labels_col', labelinfo$color)


treepdf=file.path(dbxdir, 'tree.pdf')
pdf(treepdf, h=35, w=11)
par(mar = c(3, 4, 2, 15) + 0.1)
plot(plaindend, horiz=TRUE)
dev.off()




