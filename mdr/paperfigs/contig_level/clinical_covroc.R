library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(fossil)
library(mclust)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'


####methylation distance
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))


####rank contigs by 'weirdness' of coverage
covfile=file.path(projdir, 'mdr', 'align', paste0('200708_mdr_stool16native_asmpolished.perf.sorted.cov'))
cov_cols=c('tig', 'pos', 'cov')
cov=read_tsv(covfile, col_names=cov_cols) %>%
    group_by(tig) %>%
    mutate(covfrac=cov/max(cov))
tiglens=cov %>%
    group_by(tig) %>%
    summarise(len=max(pos))
covfreq=cov %>%
    group_by(tig, cov) %>%
    summarise(freq=n()) %>%
    mutate(normfreq=freq/max(freq)) %>%
    mutate(normcov=cov/max(cov))
covpeaks=covfreq %>%
    do(findpeaks(.)) %>%
    mutate(composite=sum(h*t*s)) %>%
    arrange(-composite) %>%
    mutate(len=tiglens$len[tiglens$tig==tig]) %>%
    ungroup() %>%
    mutate(rank=dense_rank(desc(composite)))


exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
freqs=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec) %>%
    spread(key=motif, value=freq)
nacount=colSums(is.na(freqs))/dim(freqs)[1]



####add in taxonomy info
tiginfocsv=file.path(dbxdir, 'tigbins_species.csv')
tiginfo=read_csv(tiginfocsv)
chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

####roc analysis
nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])

roc=NULL
for (i in seq(0, 100, 10)) {
    if (i==0) {
        keepcov=covpeaks$tig
    } else {
        keepcov=covpeaks[-(0:i),]$tig
    }
    methchroms=methfreq %>%
        rowwise() %>%
        filter(chrom %in% keepchroms) %>%
        filter(chrom %in% keepcov)
    chrominfo=methchroms %>%
        spread(key=motif, value=freq)
    matchrominfo=as.matrix(chrominfo %>% select(-chrom))
    rownames(matchrominfo)=chrominfo$chrom
    plaindend=matchrominfo %>%
        scale %>% 
        dist %>%
        hclust %>%
        as.dendrogram
    
    truthbins=tibble(tig=labels(plaindend)) %>%
        rowwise() %>%
        filter(tig %in% tiginfo$tig) %>%
        mutate(bin=chrombins$bin[which(chrombins$rname==tig)]) %>%
        mutate(tiglen=chrombins$rlen[chrombins$rname==tig]) %>%
        filter(bin!='unknown')
    
    elimroc=get_tree_roc(plaindend, truthbins) %>%
        mutate(samp=i)

    roc=bind_rows(roc, elimroc)
}

##too crowded to include all, just do by 20s
include=seq(0, 100, 20)
rocplot=roc %>%
    filter(samp %in% include) %>%
    mutate(samp=as.character(samp)) %>%
    mutate(seqtogether=1-seqtogether) %>%
    mutate(numtogether=1-numtogether)

mycolors=colorRampPalette(brewer.pal(8, 'Set2'))(11)
covrocpdf=file.path(dbxdir, 'clinical_contig_roc_covexclude.pdf')
pdf(covrocpdf, h=8, w=11)
seqplot=ggplot(rocplot, aes(x=seqtogether, y=seqpure, colour=samp)) +
    geom_step() +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on amount of sequence') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(seqplot)
tigplot=ggplot(rocplot, aes(x=numtogether, y=numpure, colour=samp)) +
    geom_step() +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on number of contigs') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(tigplot)
dev.off()




####try length exclusion
lenroc=NULL
for (i in seq(0,50,5)) {
    methchroms=methfreq %>%
        rowwise() %>%
        filter(chrom %in% keepchroms)

    chrombinsord=chrombins %>%
        arrange(rlen) %>%
        filter(rname %in% methchroms$chrom)

    methchroms=methchroms %>%
        filter(chrom %in% chrombinsord$rname[-(1:i)])
    
    chrominfo=methchroms %>%
        spread(key=motif, value=freq)

    matchrominfo=as.matrix(chrominfo %>% select(-chrom))
    rownames(matchrominfo)=chrominfo$chrom
    plaindend=matchrominfo %>%
        scale %>% 
        dist %>%
        hclust %>%
        as.dendrogram
    
    truthbins=tibble(tig=labels(plaindend)) %>%
        rowwise() %>%
        filter(tig %in% tiginfo$tig) %>%
        mutate(bin=chrombins$bin[which(chrombins$rname==tig)]) %>%
        mutate(tiglen=chrombins$rlen[chrombins$rname==tig]) %>%
        filter(bin!='unknown')
    
    elimroc=get_tree_roc(plaindend, truthbins) %>%
        mutate(samp=i)

    lenroc=bind_rows(lenroc, elimroc)
}

lenplot=lenroc %>%
    filter(samp<=20) %>%
    mutate(samp=as.character(samp)) %>%
    mutate(seqtogether=1-seqtogether) %>%
    mutate(numtogether=1-numtogether)

covrocpdf=file.path(dbxdir, 'clinical_contig_roc_lenexclude.pdf')
pdf(covrocpdf, h=8, w=11)
seqplot=ggplot(lenplot, aes(x=seqtogether, y=seqpure, colour=samp)) +
    geom_step() +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on amount of sequence') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(seqplot)
tigplot=ggplot(lenplot, aes(x=numtogether, y=numpure, colour=samp)) +
    geom_step() +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Based on number of contigs') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(tigplot)
dev.off()
