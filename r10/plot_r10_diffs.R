library(ggplot2)
library(tidyverse)
library(gridExtra)


countfiles=list.files('/home/yfan/Dropbox/Timplab_Data/cpowgs/error', 'csv')
countsdir='/home/yfan/Dropbox/Timplab_Data/cpowgs/error/'

##plot error hists distros individually
for (i in countfiles) {
    counts=read_csv(paste0(countsdir, i), col_names=c('seq', 'freq'))
    counts=counts %>%
        arrange(-freq) %>%
        mutate(rank=c(1:dim(counts)[1]))
    
    name=substring(i, 1, nchar(i)-4)
        
    pdf(paste0(countsdir, name, '.pdf'), height=8.5, width=11)
    kmerhist=ggplot(counts, aes(x=freq)) +
        ggtitle('Error kmer histogram') +
        geom_histogram() +
        xlab('Number of errors') +
        theme_bw()
    rarefaction=ggplot(counts, aes(x=rank, y=freq)) +
        ggtitle('Error rarefaction') +
        geom_point() +
        xlab('Rank') +
        ylab('Frequency')+
        theme_bw()
    
    print(kmerhist)
    print(rarefaction)
    print(grid.table(counts[1:20,1:2], rows=NULL))
    dev.off()
}



##plot error hists and distros together
name=substring(countfiles[1], 1, nchar(countfiles[1])-4)
corr=substring(countfiles[1], 1, nchar(countfiles[1])-10)
mer=substring(countfiles[1],  nchar(countfiles[1])-8, nchar(countfiles[1])-4)

allcounts=read_csv(paste0(countsdir, countfiles[1]), col_names=c('seq', 'freq'))
allcounts=allcounts %>%
    arrange(-freq) %>%
    mutate(rank=c(1:dim(allcounts)[1])) %>%
    mutate(samp=name) %>%
    mutate(corr=corr) %>%
    mutate(mer=mer)

for (i in countfiles[2:length(countfiles)]) {
    name=substring(i, 1, nchar(i)-4)
    corr=substring(i, 1, nchar(i)-10)
    mer=substring(i, nchar(i)-8, nchar(i)-4)
    
    counts=read_csv(paste0(countsdir, i), col_names=c('seq', 'freq'))
    counts=counts %>%
        arrange(-freq) %>%
        mutate(rank=c(1:dim(counts)[1])) %>%
        mutate(samp=name) %>%
        mutate(corr=corr) %>%
        mutate(mer=mer)

    allcounts=bind_rows(allcounts, counts)
}

pdf(paste0(countsdir, 'error_profile.pdf'), height=8.5, width=11)
for (i in unique(allcounts$mer)) {
    rarefaction=ggplot(allcounts[allcounts$mer==i,], aes(x=rank, y=freq)) +
        ggtitle(i) +
        geom_point(aes(colour=corr)) +
        xlab('Rank') +
        ylab('Frequency') +
        theme_bw()
    print(rarefaction)
}
dev.off()


##tables comparing kmers across corrections
for (i in unique(allcounts$mer)) {
    mercounts=allcounts[allcounts$mer==i,] %>%
        mutate(samp=NULL) %>%
        mutate(rank=NULL)
    
    corrcomp=spread(mercounts, corr, freq) %>%
        arrange(-raw) %>%
        select(seq, raw, polished, mpolish, pilon)
    pdf(paste0(countsdir, i, '_tables.pdf'), width=5, height=10)
    print(grid.table(corrcomp[1:20,], rows=NULL))
    dev.off()

    improve=corrcomp %>%
        mutate(delta=mpolish-pilon) %>%
        arrange(-delta)
    pdf(paste0(countsdir, i, '_improvetables.pdf'), width=5, height=10)
    print(grid.table(improve[1:20,], rows=NULL))
    dev.off()
}
