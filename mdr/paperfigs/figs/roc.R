library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/figs'


####methylation distance
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))
##methfreq=readRDS(file.path(datadir, 'clin_methfreq3.rds'))
exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
methfreq=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec)

freqs=methfreq %>%
    spread(key=motif, value=freq)
nacount=colSums(is.na(freqs))/dim(freqs)[1]


####chroms to bins
chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

####add in taxonomy info
tiginfocsv=file.path(dbxdir, 'tigbins_species.csv')
tiginfo=read_csv(tiginfocsv)


####roc stuff
nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])
methchroms=methfreq %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)
chrominfo=methchroms %>%
    spread(key=motif, value=freq)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom
plaindend=matchrominfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram

##truth
truthbins=tibble(tig=labels(plaindend)) %>%
    rowwise() %>%
    filter(tig %in% tiginfo$tig) %>%
    mutate(bin=chrombins$bin[which(chrombins$rname==tig)]) %>%
    mutate(tiglen=chrombins$rlen[chrombins$rname==tig]) %>%
    filter(bin!='unknown')


library(fossil)
library(mclust)
maxheight=attributes(plaindend)$height
contams=NULL
for (height in seq(0, maxheight, .01)){
    clusts=cutree(plaindend, h=height)
    truthbins=truthbins %>%
        rowwise() %>%
        mutate(clustbin=clusts[tig])
    numclusts=length(table(truthbins$clustbin))
    
    ##togetherness - what percentage of contigs is not in the right cluster?
    ##assign bins to clusters by majority
    majclusts=truthbins %>%
        group_by(bin, clustbin) %>%
        summarise(seq=sum(tiglen)) %>%
        filter(seq==max(seq))
    ##for each contig, see if its part of the 'right' cluster
    contam=truthbins %>%
        rowwise() %>%
        mutate(inmajority=clustbin==majclusts$clustbin[majclusts$bin==bin])
    majority=sum(contam$tiglen[contam$inmajority])/sum(contam$tiglen)
    nummajority=sum(contam$inmajority)/dim(contam)[1]
    
    ##purity - what percentage of a meth cluster is not the majority bin of that cluster?
    ##assign clusters to bins according to majority composition
    pureclusts=truthbins %>%
        group_by(clustbin, bin) %>%
        summarise(seq=sum(tiglen)) %>%
        ungroup() %>%
        group_by(clustbin) %>%
        filter(seq==max(seq))
    ##for each contig, see if it's part of the 'right' cluster
    purity=truthbins %>%
        rowwise() %>%
        mutate(purity=bin==pureclusts$bin[pureclusts$clustbin==clustbin])
    pure=sum(purity$tiglen[purity$purity])/sum(purity$tiglen)
    numpure=sum(purity$purity)/dim(purity)[1]

    heightcontam=tibble(height=height,
                      numclusts=numclusts,
                      majority=majority,
                      revmajority=1-majority,
                      nummajority=nummajority,
                      revnummajority=1-nummajority,
                      pure=pure,
                      numpure=numpure)
    contams=bind_rows(contams, heightcontam)
}

ucontams=contams %>%
    select(-height) %>%
    unique() %>%
    arrange(numclusts)
percent=ucontams %>%
    select(-c(nummajority, numpure, revmajority, revnummajority)) %>%
    gather('key', 'value', -numclusts)
num=ucontams %>%
    select(-c(majority, pure, revmajority, revnummajority)) %>%
    gather('key', 'value', -numclusts)


ratings=contams %>%
    mutate(seqrate=majority*pure) %>%
    mutate(tigrate=nummajority*numpure) %>%
    unique() %>%
    select(-c(height))

allroc=NULL
for (i in 1:50) {
    randchrominfo=as.matrix(chrominfo %>% select(-chrom))
    rownames(randchrominfo)=sample(chrominfo$chrom)
    randodend=randchrominfo %>%
        scale %>% 
        dist %>%
        hclust %>%
        as.dendrogram
    randroc=get_tree_roc(randodend, truthbins) %>%
        mutate(samp=i)
    allroc=bind_rows(allroc, randroc)
}

numtips=length(chrominfo$chrom)
allrands=NULL
for (i in 1:50) {
    randchrominfo=matrix(sample(matchrominfo), nrow=dim(matchrominfo)[1], ncol=dim(matchrominfo)[2])
    rownames(randchrominfo)=chrominfo$chrom
    colnames(randchrominfo)=colnames(matchrominfo)
    randodend=randchrominfo %>%
        scale %>% 
        dist %>%
        hclust %>%
        as.dendrogram
    randroc=get_tree_roc(randodend, truthbins) %>%
        mutate(samp=i)
    allrands=bind_rows(allrands, randroc)
}


realroc=get_tree_roc(plaindend, truthbins) %>%
    mutate(samp=0)


plotallroc=bind_rows(realroc, allroc) %>%
    mutate(label=case_when(samp==0 ~ 'real', TRUE ~ 'rando')) %>%
    mutate(seqtogether=1-seqtogether) %>%
    mutate(numtogether=1-numtogether) %>%
    filter(samp<=20)
plotallrands=bind_rows(realroc, allrands) %>%
    mutate(label=case_when(samp==0 ~ 'real', TRUE ~ 'rando')) %>%
    mutate(seqtogether=1-seqtogether) %>%
    mutate(numtogether=1-numtogether) %>%
    filter(samp<=20)


plotallroccsv=file.path(dbxdir, 'plotallroc.csv')
write_csv(plotallroc, plotallroccsv)
plotallrandscsv=file.path(dbxdir, 'plotallrands.csv')
write_csv(plotallrands, plotallrandscsv)


rocpdf=file.path(dbxdir, 'roc.pdf')
pdf(randopdf, h=8, w=11)
plot1=ggplot(plotallroc %>% filter(label=='real'), aes(x=seqtogether, y=seqpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallroc %>% filter(label=='rando'), mapping=aes(x=seqtogether, y=seqpure, colour=label, alpha=.02)) +
    ggtitle('Percent Sequence: random labels, same tree structure') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot1)
plot2=ggplot(plotallrands %>% filter(label=='real'), aes(x=seqtogether, y=seqpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallrands %>% filter(label=='rando'), mapping=aes(x=seqtogether, y=seqpure, colour=label, alpha=.02)) +
    ggtitle('Percent Sequence: random structure') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot2)
plot3=ggplot(plotallroc %>% filter(label=='real'), aes(x=numtogether, y=numpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallroc %>% filter(label=='rando'), mapping=aes(x=numtogether, y=numpure, colour=label, alpha=.02)) +
    ggtitle('Percent Contigs: random labels, same tree structure') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot3)
plot4=ggplot(plotallrands %>% filter(label=='real'), aes(x=numtogether, y=numpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallrands %>% filter(label=='rando'), mapping=aes(x=numtogether, y=numpure, colour=label, alpha=.02)) +
    ggtitle('Percent Contigs: random sturcture') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot4)
dev.off()




