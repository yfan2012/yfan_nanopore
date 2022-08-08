library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(cluster)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'


invec=c('CAGAG','CCWGG', 'CMTCGAKG','CTCCAG', 'CTKVAG', 'GATC', 'GCGC', 'GCWGC', 'GGCC', 'GGNNCC', 'GGWCC', 'TCCGGA', 'GCCGGC', 'RGCGCY')
methfreq=readRDS(file.path(datadir, 'methfreq.rds')) %>%
    rowwise() %>%
    filter(motif %in% invec)

chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])
methchroms=methfreq %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)
chrominfo=methchroms %>%
    spread(key=motif, value=freq) %>%
    filter(chrom %in% chrombins$rname)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom
plaindend=matchrominfo %>%
    scale %>% 
    dist %>%
    diana %>%
    as.dendrogram

truthbins=tibble(tig=labels(plaindend)) %>%
    rowwise() %>%
    ##filter(tig %in% tiginfo$tig) %>%
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

library(foreach)
library(doParallel)
registerDoParallel(25)

allroc=foreach(i=1:50, .combine=rbind) %dopar% {
    cat(paste0(as.character(i),"\n"), file=stdout())
    randchrominfo=as.matrix(chrominfo %>% select(-chrom))
    rownames(randchrominfo)=sample(chrominfo$chrom)
    randodend=randchrominfo %>%
        scale %>% 
        dist %>%
        diana %>%
        as.dendrogram
    randroc=get_tree_roc(randodend, truthbins) %>%
        mutate(samp=i)
    return(randroc)
}


numtips=length(chrominfo$chrom)
allrands=foreach(i=1:50, .combine=rbind) %dopar% {
    cat(paste0(as.character(i),"\n"), file=stdout())
    randchrominfo=matrix(sample(matchrominfo), nrow=dim(matchrominfo)[1], ncol=dim(matchrominfo)[2])
    rownames(randchrominfo)=chrominfo$chrom
    colnames(randchrominfo)=colnames(matchrominfo)
    randodend=randchrominfo %>%
        scale %>% 
        dist %>%
        diana %>%
        as.dendrogram
    randroc=get_tree_roc(randodend, truthbins) %>%
        mutate(samp=i)
    return(randroc)
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

plotallroccsv=file.path(dbxdir, 'plotallroc_diana.csv')
write_csv(plotallroc, plotallroccsv)
plotallrandscsv=file.path(dbxdir, 'plotallrands_diana.csv')
write_csv(plotallrands, plotallrandscsv)


plotallroc=read_csv(file.path(dbxdir, 'plotallroc_diana.csv')) %>%
    mutate(randtype='leaf') %>%
    filter(samp<=5)
plotallrands=read_csv(file.path(dbxdir, 'plotallrands_diana.csv')) %>%
    mutate(randtype='tree') %>%
    filter(samp<=5) %>%
    filter(label=='rando')
rocplot=bind_rows(plotallroc, plotallrands) %>%
    mutate(samp=as.character(samp))


rocpdf=file.path(dbxdir, 'roc_diana.pdf')
pdf(rocpdf, h=8, w=11)
seqplot=ggplot(rocplot %>% filter(label=='real'), aes(x=seqtogether, y=seqpure, colour=samp)) +
    geom_step() +
    geom_point(rocplot %>% filter(label=='rando'), mapping=aes(x=seqtogether, y=seqpure, colour=randtype, alpha=.02)) +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on amount of sequence') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(seqplot)
tigplot=ggplot(rocplot %>% filter(label=='real'), aes(x=numtogether, y=numpure, colour=samp)) +
    geom_step() +
    geom_point(rocplot %>% filter(label=='rando'), mapping=aes(x=numtogether, y=numpure, colour=randtype, alpha=.02)) +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on number of contigs') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(tigplot)
dev.off()
