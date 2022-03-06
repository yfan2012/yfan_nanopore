library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'barnyard/contig_level')
dbxdir='~/gdrive/mdr/paperfigs/contig_level'


get_mat <- function (sampname, picked) {
    methfile=file.path(datadir, paste0(sampname, '.barocdes_methcalls.csv'))
    methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
    
    meth=read_csv(methfile, col_names=methcols) %>%
        group_by(chrom, pos, strand, motif) %>%
        summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
        mutate(methfrac=methnum/(methnum+umethnum))
    
    cluster_copy(cluster, 'findMethFreq')
    
    methgrouped=meth %>%
        filter(sum(methnum+umethnum)>5) %>%
        ##filter(sum(methnum+umethnum)>1) %>%
        group_by(chrom, motif) %>%
        partition(cluster)
    methfreq=methgrouped %>%
        do(findMethFreq(.))  %>%
        collect() %>%
        summarise(freq=mean(methfrac))
    
    nummotifs=length(table(methfreq$motif))
    if (length(picked) >1){
        keepmotifs=picked
    } else {
        keepmotifs=names(table(methfreq$motif)[table(methfreq$motif)==3])
    }
    ##exclude based on clinical analysis
    exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
    
    submethfreq=methfreq %>%
        filter(motif %in% keepmotifs) %>%
        filter(!motif %in% exvec) %>%
        filter(motif != 'GCGC')
    
    return(submethfreq)
}

get_dists <- function(submethfreq) {
    freqs=submethfreq %>%
        spread(key=motif, value=freq) %>%
        ungroup()
    
    ##make heatmap
    matfreqs=as.matrix(freqs %>% select(-chrom))
    rownames(matfreqs)=freqs$chrom
    freqdists=as.matrix(dist(matfreqs))
        
    freqdistsdf=as_tibble(freqdists) %>%
        mutate(chroms=rownames(freqdists)) %>%
        gather(key=chroms2, value=dists, -chroms) %>%
        mutate(rounded=round(dists,2))
    return(freqdistsdf)
}

prefix='211005_mdr_barnyard_mix'
samp=c('5','6')
sampname=paste0(prefix, samp)

distinfo=NULL
blind=c()
for (i in sampname) {
    sampdistinfo=get_dists(i, blind) %>%
        mutate(samp=i)
    distinfo=bind_rows(distinfo, sampdistinfo)
}

pickinfo=NULL
picked=c('CCWGG', 'GATC', 'KCCGGM')
for (i in sampname) {
    samppickinfo=get_dists(i, picked) %>%
        mutate(samp=i)
    pickinfo=bind_rows(pickinfo, samppickinfo)
}
    

##plot heatmap
distplotspdf=file.path(dbxdir, paste0('barnyard_heatmaps.pdf'))
pdf(distplotspdf, h=6, w=9)
for (i in sampname) {
    sampdistinfo=distinfo %>%
        filter(samp==i)
    plot=ggplot(sampdistinfo, aes(x=chroms, y=chroms2)) +
        geom_tile(aes(fill = dists)) +
        geom_text(aes(label = rounded)) +
        scale_fill_gradient(low = "white", high = "red") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(plot)
}
dev.off()
    

pickplotspdf=file.path(dbxdir, paste0('barnyard_heatmaps_picked.pdf'))
pdf(pickplotspdf, h=6, w=9)
for (i in sampname) {
    samppickinfo=pickinfo %>%
        filter(samp==i)
    plot=ggplot(samppickinfo, aes(x=chroms, y=chroms2)) +
        geom_tile(aes(fill = dists)) +
        geom_text(aes(label = rounded)) +
        scale_fill_gradient(low = "white", high = "red") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(plot)
}
dev.off()
