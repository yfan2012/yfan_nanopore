library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_asm'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'


####plasmids - figure out which contigs need to be included
binplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.hiC.plasmidfinder.tsv')
tigplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.plasmidfinder.tsv')
plas_cols=c('file', 'seq', 'start', 'end', 'strand', 'gene', 'coverage', 'covmap', 'gaps', 'covfrac', 'ident', 'db', 'acc', 'prod', 'res')
binplas=read_tsv(binplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -coverage, -covmap, -covmap, -db, -prod, -res) %>%
    filter(ident>=95)
tigplas=read_tsv(tigplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -coverage, -covmap, -covmap, -db, -prod, -res) %>%
    filter(ident>=95) ##%>%
    ##rowwise() %>%
    ##filter(gene %in% binplas$gene) %>%
    ##mutate(bin=binplas$seq[binplas$gene==gene])
plastigs=names(table(tigplas$seq))
missingplas=binplas$gene[!binplas$gene %in% tigplas$gene]




####read meth data
methfile=file.path(datadir, 'clin_barocdes_methcalls.csv')
methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
meth=read_csv(methfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))




####try no restricitons on just the dropped plasmid sequences
cluster_copy(cluster, 'findMethFreq')
methgrouped=meth %>%
    filter(sum(methnum+umethnum)>5) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methloci=methgrouped %>%
    do(findMethFreq(.))  %>%
    collect()
methfreq=methloci %>%
    summarise(freq=mean(methfrac))

plasgrouped=meth %>%
    filter(chrom %in% plastigs) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
plasloci=plasgrouped %>%
    do(findMethFreq(.))  %>%
    collect()
plasfreq=plasloci %>%
    summarise(freq=mean(methfrac))


allmotifs=names(table(methfreq$motif))

##check how incomplete some of the plasmid contigs are
##looks like most important plasmids need to have 
tigcomplete=methfreq %>%
    summarise(count=n())

plastiginfo=tigplas %>%
    rowwise() %>%
    mutate(count=tigcomplete$count[tigcomplete$chrom==seq])


imputezero <- function(test, mindata) {
    ##fill in missing values with zeros
    ##require that mindata points are already filled
    if (dim(test)[1]>=mindata) {
        missingmotif=allmotifs[!allmotifs %in% test$motif]
        missinginfo=tibble(chrom=test$chrom[1], motif=missingmotif, freq=0)
        fullinfo=bind_rows(test, missinginfo)
    } else {
        fullinfo=test
    }
    return(fullinfo)
}

methchroms=methfreq %>%
    do(imputezero(., 9))

nummotifs=length(allmotifs)
keepchroms=names(table(methchroms$chrom)[table(methchroms$chrom)==nummotifs])
methchroms=methchroms %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)
    
##make matrix of meth info by contig
chrominfo=methchroms %>%
    spread(key=motif, value=freq)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom

##make dstance matrix
chromdists=as_tibble(as.matrix(dist(matchrominfo))) %>%
    mutate(chroms=chrominfo$chrom) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist, 2))

####add in bin info from mummer
chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)
for (i in names(table(methchroms$chrom))) {
    if (!(i %in% chrombins$rname)) {
        info=tibble(rname=i, bin='unknown', total_rcov=NA)
        chrombins=bind_rows(chrombins, info)
    }
}
methbins=methchroms %>%
    rowwise() %>%
    mutate(bin=chrombins$bin[chrombins$rname==chrom])

##set up colors
binnames=c()
for (i in rownames(matchrominfo)) {
    bin=chrombins$bin[chrombins$rname==i]
    binnames=c(binnames, bin)
}
numcolors=length(names(table(binnames)))-1
mycolors=c(colorRampPalette(brewer.pal(8, 'Set2'))(numcolors), '#000000')
colorkey=tibble(contig=rownames(table(binnames)), color=mycolors)
bincolors=c()
for (i in binnames) {
    color=colorkey$color[colorkey$contig==i]
    bincolors=c(bincolors, color)
}
bincolorinfo=data.frame(color=bincolors)
rownames(bincolorinfo)=rownames(matchrominfo)

dend=matchrominfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
label_order=labels(dend)
label_colors=bincolorinfo[label_order,]
dend=dend %>%
    set('labels_col', label_colors)

chrombinspdf=file.path(dbxdir, 'clinical_contig_clusters_imputezero.pdf')
pdf(chrombinspdf, h=40, w=8)
par(mar = c(5, 4, 2, 10) + 0.1)
plot(dend, horiz=TRUE)
dev.off()

binnedinfo=matchrominfo[binnames!='unknown',]
binneddend=binnedinfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
label_order=labels(binneddend)
label_colors=bincolorinfo[label_order,]
binneddend=binneddend %>%
    set('labels_col', label_colors)

chrombinspdf=file.path(dbxdir, 'clinical_contig_clusters_imputezero_binned.pdf')
pdf(chrombinspdf, h=45, w=11)
par(mar = c(5, 4, 2, 10) + 0.1)
plot(binneddend, horiz=TRUE)
dev.off()




