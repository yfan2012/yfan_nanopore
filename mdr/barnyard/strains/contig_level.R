library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'barnyard/strains')
prefix='220131_mdr_barnyard_'

dbxdir='~/gdrive/mdr/paperfigs/contig_level'


methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
staphmethfile=file.path(datadir, 'contig_level', paste0(prefix, 'st3294.barocdes_methcalls.csv'))
ecolimethfile=file.path(datadir, 'contig_level', paste0(prefix, 'st3689.barocdes_methcalls.csv'))
staphmeth=read_csv(staphmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
ecolimeth=read_csv(ecolimethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))


cluster_copy(cluster, 'findMethFreq')
getFreqdf <- function(meth) {
    methgrouped=meth %>%
        filter(sum(methnum+umethnum)>5) %>%
        ##filter(sum(methnum+umethnum)>1) %>%
        group_by(chrom, motif) %>%
        partition(cluster)
    methfreq=methgrouped %>%
        do(findMethFreq(.))  %>%
        collect() %>%
        summarise(freq=mean(methfrac))
}
staphfreq=getFreqdf(staphmeth)
ecolifreq=getFreqdf(ecolimeth)
freq=bind_rows(staphfreq %>% mutate(samp='staph'), ecolifreq %>% mutate(samp='ecoli'))
keepchroms=names(table(freq$motif))[table(freq$motif)==4]
##get rid of gcgc because it's sort of redundant with rgcgcy
freq=freq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC')


splasfreq=freq %>%
    filter(!(chrom=='PRW62' & samp=='ecoli'))
eplasfreq=freq %>%
    filter(!(chrom=='PRW62' & samp=='staph'))


getFreqMatrix <- function(freq) {
    chrominfo=freq %>%
        spread(key=motif, value=freq) %>%
        ungroup()
    matchrominfo=as.matrix(chrominfo %>% select(-chrom, -samp))
    rownames(matchrominfo)=paste0(chrominfo$samp, '_', chrominfo$chrom)
    return(matchrominfo)
}

freqmat=getFreqMatrix(freq)
freqmat=rbind(freqmat, (freqmat['staph_PRW62',]+freqmat['ecoli_PRW62',])/2)
rownames(freqmat)[5]='both_PRW62'
splasmat=getFreqMatrix(splasfreq)
eplasmat=getFreqMatrix(eplasfreq)

getDists <- function(freqmat) {
    dists=as_tibble(as.matrix(dist(freqmat))) %>%
        mutate(chroms=rownames(freqmat)) %>%
        gather(key=chroms2, value=dist, -chroms) %>%
        mutate(rounded=round(dist, 2))
    return(dists)
}

dists=getDists(freqmat)
splasdists=getDists(splasmat)
eplasdists=getDists(eplasmat)


distpdf=file.path(dbxdir, 'barnyard_strain_distances.pdf')
pdf(distpdf, h=6, w=6)
plot=ggplot(dists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
dev.off()




samps=c('both_plas', 'staph_plas', 'ecoli_plas')
covcols=c('chrom', 'pos', 'cov')
cov=NULL
for (i in samps) {
    covfile=file.path(datadir, 'align', paste0(i, '.genomecov'))
    sampcov=read_tsv(covfile, covcols) %>%
        group_by(chrom) %>%
        summarise(meancov=mean(cov)) %>%
        mutate(samp=i)
    cov=bind_rows(cov, sampcov)
}
runs=c('220131_mdr_barnyard_st3294', '220131_mdr_barnyard_st3689')
runcov=NULL
for (i in runs) {
    covfile=file.path(datadir, 'align', paste0(i, '.genomecov'))
    sampcov=read_tsv(covfile, covcols) %>%
        group_by(chrom) %>%
        summarise(meancov=mean(cov)) %>%
        mutate(samp=i)
    runcov=bind_rows(runcov, sampcov)
}


splasmethfile=file.path(datadir, 'contig_level', 'staph_plas.barocdes_methcalls.csv')
eplasmethfile=file.path(datadir, 'contig_level', 'ecoli_plas.barocdes_methcalls.csv')
aplasmethfile=file.path(datadir, 'contig_level', 'both_plas.barocdes_methcalls.csv')
splasmeth=read_csv(splasmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
eplasmeth=read_csv(eplasmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
aplasmeth=read_csv(aplasmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))


splasfreq=getFreqdf(splasmeth)
eplasfreq=getFreqdf(eplasmeth)
aplasfreq=getFreqdf(aplasmeth)

##get rid of gcgc because it's sort of redundant with rgcgcy
sfreq=splasfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    mutate(samp='s')
efreq=eplasfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    mutate(samp='e')
afreq=aplasfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    mutate(samp='a')

sfreqmat=getFreqMatrix(sfreq)
rownames(sfreqmat)=c('chr_ecoli', 'chr_staph', 'plas_staph')
efreqmat=getFreqMatrix(efreq)
rownames(efreqmat)=c('chr_ecoli', 'chr_staph', 'plas_ecoli')
afreqmat=getFreqMatrix(afreq)
rownames(afreqmat)=c('chr_ecoli', 'chr_staph', 'plas_both')

sdists=getDists(sfreqmat)
edists=getDists(efreqmat)
adists=getDists(afreqmat)

combine_freqmat=rbind(sfreqmat, efreqmat[3,], afreqmat[3,])
rownames(combine_freqmat)=c('ecoli_chr', 'staph_chr', 'staph_plas', 'ecoli_plas', 'all_plas')

cdists=getDists(combine_freqmat)

sepdistpdf=file.path(dbxdir, 'barnyard_strain_distances_sep.pdf')
pdf(sepdistpdf, h=6, w=6)
splot=ggplot(sdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(splot)
eplot=ggplot(edists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(eplot)
aplot=ggplot(adists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(aplot)
cplot=ggplot(cdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(cplot)
dev.off()










freq=bind_rows(staphfreq %>% mutate(samp='staph'), ecolifreq %>% mutate(samp='ecoli'))
keepchroms=names(table(freq$motif))[table(freq$motif)==4]

freq=freq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC')


splasfreq=freq %>%
    filter(!(chrom=='PRW62' & samp=='ecoli'))
eplasfreq=freq %>%
    filter(!(chrom=='PRW62' & samp=='staph'))
