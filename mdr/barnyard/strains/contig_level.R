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
        summarise(freq=mean(methfrac), numloci=n())
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


readcov=NULL
for (i in samps) {
    covfile=file.path(datadir, 'align', paste0(i, '.genomecov'))
    sampcov=read_tsv(covfile, covcols) %>%
        mutate(samp=i)
    readcov=bind_rows(readcov, sampcov)
}


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


splasfreq=getFreqdf(splasmeth)
eplasfreq=getFreqdf(eplasmeth)
aplasfreq=getFreqdf(aplasmeth)


##get rid of gcgc because it's sort of redundant with rgcgcy
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

abbrevpdf=file.path(dbxdir, 'barnyard_strain_fewermotifs.pdf')
pdf(abbrevpdf, h=6, w=6)
aplot=ggplot(tdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label=rounded)) +
    scale_fill_gradient(low='white', high='red') +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1))
print(aplot)
dev.off()

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




####figure out if there's a good coverage based way to exclude some motifs
smethcov=full_join(splasmeth, readcov %>% filter(samp=='staph_plas')) %>%
    na.omit %>%
    mutate(methcov=(methnum+umethnum)/cov) %>%
    mutate(methtot=methnum/cov)
emethcov=full_join(eplasmeth, readcov %>% filter(samp=='ecoli_plas')) %>%
    na.omit %>%
    mutate(methcov=(methnum+umethnum)/cov) %>%
    mutate(methtot=methnum/cov)
amethcov=full_join(aplasmeth, readcov %>% filter(samp=='both_plas')) %>%
    na.omit %>%
    mutate(methcov=(methnum+umethnum)/cov) %>%
    mutate(methtot=methnum/cov)

methcov=bind_rows(smethcov, emethcov, amethcov)

methcovpdf=file.path(dbxdir, 'barnyard_strain_methcov.pdf')
pdf(methcovpdf, w=11, h=8)
dev.off()



####test totcall stuff
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













splasonly=splasloci %>% filter(chrom=='PRW62') %>%
    mutate(org='staph')
eplasonly=eplasloci %>% filter(chrom=='PRW62') %>% 
    mutate(org='ecoli')
plasonly=bind_rows(splasonly, eplasonly)

plascomp=eplasonly
plascomp$staph=splasonly$totcalls
    

tshootpdf=file.path(dbxdir, 'barnyard_troubleshoot.pdf')
pdf(tshootpdf, w=11, h=8)
plot=ggplot(plasonly, aes(x=totcalls, colour=org, fill=org, alpha=.2)) +
    geom_density() +
    facet_wrap(~motif) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
print(plot)
scatter=ggplot(plascomp, aes(x=totcalls, y=staph, alpha=.1)) +
    geom_point() +
    facet_wrap(~motif) +
    xlim(0,250) +
    ylim(0,250) +
    geom_abline(slope=1, intercept=0) +
    theme_bw()
print(scatter)
dev.off()


source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')
cluster_copy(cluster, 'findMethFreqtest')
testfreqdf <- function(meth) {
    methgrouped=meth %>%
        ##filter(sum(methnum+umethnum)>5) %>%
        ##filter(sum(methnum+umethnum)>1) %>%
        group_by(chrom, motif) %>%
        partition(cluster)
    methfreq=methgrouped %>%
        do(findMethFreqtest(.))  %>%
        collect() %>%
        summarise(freq=mean(methfrac), numloci=n())
}

stestfreq=testfreqdf(splasmeth)
etestfreq=testfreqdf(eplasmeth)
atestfreq=testfreqdf(aplasmeth)

cluster_copy(cluster, 'findMethFreqcov')
cluster_copy(cluster, 'readcov')
collapse_cov <- function(meth, cov) {
    methgrouped=meth %>%
        ##filter(sum(methnum+umethnum)>5) %>%
        ##filter(sum(methnum+umethnum)>1) %>%
        group_by(chrom, motif) %>%
        partition(cluster)
    methfreq=methgrouped %>%
        do(findMethFreqcov(., cov))  %>%
        collect() ##%>%
        ##summarise(freq=mean(methfrac), numloci=n())
}


stestfreq=collapse_cov(splasmeth, readcov %>% filter(samp=='staph_plas'))
etestfreq=collapse_cov(eplasmeth, readcov %>% filter(samp=='ecoli_plas'))
atestfreq=collapse_cov(aplasmeth, readcov %>% filter(samp=='both_plas'))
saveRDS(stestfreq, 'stestfreq_cov.rds')
saveRDS(etestfreq, 'etestfreq_cov.rds')
saveRDS(atestfreq, 'atestfreq_cov.rds')


stestfreq=stestfreq %>%
    mutate(methfrac=methnum/

stest=stestfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='s')
etest=etestfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='s')
atest=atestfreq %>%
    filter(motif %in% keepchroms) %>%
    filter(motif!='GCGC') %>%
    select(-numloci) %>%
    mutate(samp='s')

stestmat=getFreqMatrix(stest)
rownames(stestmat)=c('chr_ecoli', 'chr_staph', 'plas_staph')
etestmat=getFreqMatrix(etest)
rownames(etestmat)=c('chr_ecoli', 'chr_staph', 'plas_ecoli')
atestmat=getFreqMatrix(atest)
rownames(atestmat)=c('chr_ecoli', 'chr_staph', 'plas_both')

stestdists=getDists(stestmat)
etestdists=getDists(etestmat)
atestdists=getDists(atestmat)

combine_testmat=rbind(stestmat, etestmat[3,], atestmat[3,])
rownames(combine_testmat)=c('ecoli_chr', 'staph_chr', 'staph_plas', 'ecoli_plas', 'all_plas')
ctests=getDists(combine_testmat)

fewloci=unique(stestfreq$motif[stestfreq$numloci<10]) ##same for all samps
ttestmat=combine_testmat[,-which(colnames(combine_testmat) %in% fewloci)]
ttests=getDists(ttestmat)
