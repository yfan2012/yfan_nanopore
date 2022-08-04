library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(18)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'barnyard/strains')
prefix='220131_mdr_barnyard_'

dbxdir='~/gdrive/mdr/paperfigs/figs'



cluster_copy(cluster, 'findMethFreq')
getFreqMatrix <- function(freq) {
    chrominfo=freq %>%
        spread(key=motif, value=freq) %>%
        ungroup()
    matchrominfo=as.matrix(chrominfo %>% select(-chrom, -samp))
    rownames(matchrominfo)=paste0(chrominfo$samp, '_', chrominfo$chrom)
    return(matchrominfo)
}

getDists <- function(freqmat) {
    dists=as_tibble(as.matrix(dist(freqmat))) %>%
        mutate(chroms=rownames(freqmat)) %>%
        gather(key=chroms2, value=dist, -chroms) %>%
        mutate(rounded=round(dist, 2))
    return(dists)
}

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

readcov=NULL
for (i in samps) {
    covfile=file.path(datadir, 'align', paste0(i, '.genomecov'))
    sampcov=read_tsv(covfile, covcols) %>%
        mutate(samp=i)
    readcov=bind_rows(readcov, sampcov)
}


methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
splasmethfile=file.path(datadir, 'contig_level', 'staph_plas.barocdes_methcalls.csv')
eplasmethfile=file.path(datadir, 'contig_level', 'ecoli_plas.barocdes_methcalls.csv')
aplasmethfile=file.path(datadir, 'contig_level', 'both_plas.barocdes_methcalls.csv')
scov=readcov %>% filter(samp=='staph_plas')
splasmeth=read_csv(splasmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
splasmeth=inner_join(splasmeth, scov %>% select(-samp), by=c('chrom', 'pos'))
ecov=readcov %>% filter(samp=='ecoli_plas')
eplasmeth=read_csv(eplasmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
eplasmeth=inner_join(eplasmeth, ecov %>% select(-samp), by=c('chrom', 'pos'))
acov=readcov %>% filter(samp=='both_plas')
aplasmeth=read_csv(aplasmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
aplasmeth=inner_join(aplasmeth, acov %>% select(-samp), by=c('chrom', 'pos'))


collapse_loci <- function(meth) {
    methgrouped=meth %>%
        group_by(chrom, motif) %>%
        partition(cluster)
    methfreq=methgrouped %>%
        do(findMethFreq(.))  %>%
        collect()
}

splasloci=collapse_loci(splasmeth)
eplasloci=collapse_loci(eplasmeth)
aplasloci=collapse_loci(aplasmeth)


cov_total <- function(plasloci) {
    testfreq=plasloci %>%
        ungroup() %>%
        mutate(methfrac=methnum/cov) %>%
        group_by(motif, chrom) %>%
        summarise(freq=mean(methfrac), numloci=n())
    return(testfreq)
}

splasfreq=cov_total(splasloci)
eplasfreq=cov_total(eplasloci)
aplasfreq=cov_total(aplasloci)

allplasfreq=bind_rows(splasfreq, eplasfreq, aplasfreq)
keepchroms=names(table(allplasfreq$motif)[table(allplasfreq$motif)==9])


sfreq=splasfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='s')
efreq=eplasfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='e')
afreq=aplasfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
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

barcode=c('GATC', 'CCWGG', 'CAGAG', 'CTKVAG', 'GTWWAC')

##fewloci=unique(splasfreq$motif[splasfreq$numloci<10]) ##same for all samps
fewloci=splasfreq$motif[!splasfreq$motif %in% barcode]
tfreqmat=combine_freqmat[,-which(colnames(combine_freqmat) %in% fewloci)]
tdists=getDists(tfreqmat)

covtotpdf=file.path(dbxdir, 'barnyard_heatmaps.pdf')
pdf(covtotpdf, h=6, w=6)
tplot=ggplot(tdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "#66C2A5", high = "white") +
    ggtitle('Selected motifs') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(tplot)
dev.off()



####Show non-normalized mix
staphmethfile=file.path(datadir, 'contig_level', '220131_mdr_barnyard_st3294.barocdes_methcalls.csv')
staphplasmeth=read_csv(staphmethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
ecolimethfile=file.path(datadir, 'contig_level', '220131_mdr_barnyard_st3689.barocdes_methcalls.csv')
ecoliplasmeth=read_csv(ecolimethfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))
mixplasmeth=full_join(x=ecoliplasmeth %>% filter(chrom=='PRW62'),
                      y=staphplasmeth %>% filter(chrom=='PRW62'),
                      by=c('chrom', 'pos', 'strand', 'motif')) %>%
    mutate(methnum=methnum.x+methnum.y) %>%
    mutate(umethnum=umethnum.x+umethnum.y) %>%
    select(c(chrom, pos, strand, motif, methnum, umethnum)) %>%
    mutate(methfrac=methnum/(methnum+umethnum))


staphloci=collapse_loci(staphplasmeth)
ecoliloci=collapse_loci(ecoliplasmeth)
mixloci=collapse_loci(mixplasmeth)

get_freq_tib <- function(plasloci) {
    testfreq=plasloci %>%
        ungroup() %>%
        group_by(motif, chrom) %>%
        summarise(freq=mean(methfrac), numloci=n())
    return(testfreq)
}

staphfreq=get_freq_tib(staphloci)
ecolifreq=get_freq_tib(ecoliloci)
mixfreq=get_freq_tib(mixloci)

staphfilt=staphfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='s')
ecolifilt=ecolifreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='e')
mixfilt=mixfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='a')

staphmat=getFreqMatrix(staphfilt)
ecolimat=getFreqMatrix(ecolifilt)
mixmat=getFreqMatrix(mixfilt)

freqmat=rbind(staphmat, ecolimat, mixmat)
rownames(freqmat)=c('staph_chr', 'staph_plas', 'ecoli_chr', 'ecoli_plas', 'all_plas')
filtfreqmat=freqmat[,-which(colnames(freqmat) %in% fewloci)]
filtdists=getDists(filtfreqmat)


##show non-coverage based meth calls
snocovfreq=get_freq_tib(splasloci)
enocovfreq=get_freq_tib(eplasloci)
anocovfreq=get_freq_tib(aplasloci)

snfreq=snocovfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='s')
enfreq=enocovfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='e')
anfreq=anocovfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='a')

snfreqmat=getFreqMatrix(snfreq)
rownames(sfreqmat)=c('chr_ecoli', 'chr_staph', 'plas_staph')
enfreqmat=getFreqMatrix(enfreq)
rownames(efreqmat)=c('chr_ecoli', 'chr_staph', 'plas_ecoli')
anfreqmat=getFreqMatrix(anfreq)
rownames(afreqmat)=c('chr_ecoli', 'chr_staph', 'plas_both')

sndists=getDists(snfreqmat)
endists=getDists(enfreqmat)
andists=getDists(anfreqmat)

nfreqmat=rbind(snfreqmat, enfreqmat[3,], anfreqmat[3,])
rownames(nfreqmat)=c('ecoli_chr', 'staph_chr', 'staph_plas', 'ecoli_plas', 'all_plas')
tnfreqmat=nfreqmat[,-which(colnames(combine_freqmat) %in% fewloci)]
tndists=getDists(tnfreqmat)



covtotpdf=file.path(dbxdir, 'barnyard_heatmaps_unnormalized.pdf')
pdf(covtotpdf, h=6, w=6)
tplot=ggplot(filtdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    ggtitle('No normalization') +
    scale_fill_gradient(low = "#66C2A5", high = "white") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(tplot)
tnplot=ggplot(tndists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    ggtitle('Coverage depth normalization') +
    scale_fill_gradient(low = "#66C2A5", high = "white") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(tnplot)
dev.off()

