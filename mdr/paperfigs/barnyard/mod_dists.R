library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(Biostrings)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'barnyard/strains')
prefix='220131_mdr_barnyard_'

dbxdir='~/gdrive/mdr/paperfigs/contig_level'



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


reffile=file.path(projdir, 'barnyard/ref/allrefs.fa')
ref=readDNAStringSet(reffile, format='fasta')
names(ref)=c('CP017100.1', 'PRW62', 'NC_007795.1')


dcmlocipos=as_tibble(start(vmatchPattern('CCWGG', ref, fixed=FALSE))) %>%
    select(-group)
dcmlocineg=as_tibble(start(vmatchPattern('CCWGG', ref, fixed=FALSE))) %>%
    select(-group) %>%
    mutate(value=value+2)
dcmloci=bind_rows(dcmlocipos, dcmlocineg)
colnames(dcmloci)=c('chrom', 'pos')
damlocipos=as_tibble(start(vmatchPattern('GATC', ref, fixed=FALSE))) %>%
    select(-group)
damlocineg=as_tibble(start(vmatchPattern('GATC', ref, fixed=FALSE))) %>%
    mutate(value=value+1) %>%
    select(-group)
damloci=bind_rows(damlocipos, damlocineg)
colnames(damloci)=c('chrom', 'pos')


splasdcm=left_join(dcmloci, splasmeth) %>%
    mutate(samp='staph') %>%
    filter(chrom!='CP017100.1')
eplasdcm=left_join(dcmloci, eplasmeth) %>%
    mutate(samp='ecoli') %>%
    filter(chrom!='NC_007795.1')
aplasdcm=left_join(dcmloci, aplasmeth) %>%
    mutate(samp='both') %>%
    filter(chrom=='PRW62')
dcm=bind_rows(splasdcm, eplasdcm, aplasdcm) %>%
    filter(motif=='CCWGG') %>%
    mutate(methcall=methnum/cov) %>%
    mutate(umethcall=umethnum/cov) %>%
    mutate(unkcall=(cov-methnum-umethnum)/cov) %>%
    select(-c(pos, strand, motif, methfrac, cov, methnum, umethnum)) %>%
    gather(key=call, value=value, -samp, -chrom) %>%
    mutate(name=paste0(chrom, '_', samp))


splasdam=left_join(damloci, splasmeth) %>%
    mutate(samp='staph') %>%
    filter(chrom!='CP017100.1')
eplasdam=left_join(damloci, eplasmeth) %>%
    mutate(samp='ecoli') %>%
    filter(chrom!='NC_007795.1')
aplasdam=left_join(damloci, aplasmeth) %>%
    mutate(samp='both') %>%
    filter(chrom=='PRW62')
dam=bind_rows(splasdam, eplasdam, aplasdam) %>%
    filter(motif=='GATC') %>%
    mutate(methcall=methnum/cov) %>%
    mutate(umethcall=umethnum/cov) %>%
    mutate(unkcall=(cov-methnum-umethnum)/cov) %>%
    select(-c(pos, strand, motif, methfrac, cov, methnum, umethnum)) %>%
    gather(key=call, value=value, -samp, -chrom) %>%
    mutate(name=paste0(chrom, '_', samp))


distpdf=file.path(dbxdir, 'methdists.pdf')
pdf(distpdf, h=5, w=11)
damplot=ggplot(dam, aes(y=value, x=chrom, color=name, fill=name, alpha=.1)) +
    geom_violin(position='identity') +
    facet_wrap(~call, scales='free') +
    ggtitle('GATC') +
    theme_bw()
dcmplot=ggplot(dcm, aes(y=value, x=chrom, color=name, fill=name, alpha=.1)) +
    geom_violin(position='identity') +
    facet_wrap(~call, scales='free') +
    ggtitle('CCWGG') +
    theme_bw()
print(damplot)
print(dcmplot)
dev.off()
