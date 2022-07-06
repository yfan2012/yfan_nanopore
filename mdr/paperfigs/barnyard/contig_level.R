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

dbxdir='~/gdrive/mdr/paperfigs/contig_level'



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
        summarise(freq=mean(methfrac), numloci=n())
}

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
runcovcsv='~/gdrive/mdr/paperfigs/figs/barnyard_runcov.csv'
write_csv(runcov, runcovcsv)

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
        ##filter(sum(methnum+umethnum)>5) %>%
        ##filter(sum(methnum+umethnum)>1) %>%
        group_by(chrom, motif) %>%
        partition(cluster)
    methfreq=methgrouped %>%
        do(findMethFreq(.))  %>%
        collect() ##%>%
        ##summarise(freq=mean(methfrac), numloci=n())
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
cdists=getDists(combine_freqmat)

fewloci=unique(splasfreq$motif[splasfreq$numloci<10]) ##same for all samps
tfreqmat=combine_freqmat[,-which(colnames(combine_freqmat) %in% fewloci)]
tdists=getDists(tfreqmat)


covtotpdf=file.path(dbxdir, 'barnyard_strain_dists_cov.pdf')
pdf(covtotpdf, h=6, w=6)
cplot=ggplot(cdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "#FC8D62", high = "white") +
    ggtitle('All motifs') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(cplot)
cplot=ggplot(cdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "#FC8D62") +
    ggtitle('All motifs') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(cplot)
tplot=ggplot(tdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "#FC8D62", high = "white") +
    ggtitle('Selected motifs') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(tplot)
tplot=ggplot(tdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = 'white', high = "#FC8D62") +
    ggtitle('Selected motifs') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(tplot)
dev.off()




####unmeth check
chrmeth=splasmeth %>%
    filter(chrom!='PRW62') %>%
    mutate(org=case_when(chrom=='CP017100.1' ~ 'ecoli',
                         TRUE ~ 'staph'))
methscatterpdf=file.path(dbxdir, 'barnyard_meth_scatter.pdf')
pdf(methscatterpdf, h=8, w=15)
for (i in c('CCWGG', 'GATC')) {
    plotchr=chrmeth %>%
        filter(motif==i)
    plot=ggplot(plotchr, aes(x=umethnum, y=methnum, colour=org, alpha=.05)) +
        geom_point() +
        facet_wrap(~org) +
        ggtitle(i) +
        scale_colour_brewer(palette='Set2') +
        theme_bw()
    print(plot)
}
dev.off()
