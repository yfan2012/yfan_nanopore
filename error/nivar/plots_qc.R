library(ggplot2)
library(tidyverse)
library(tidyr)
library(gridExtra)
library(doParallel)
library(R.utils)

dbxdir='~/Dropbox/yfan/nivar/'

##plot busco
buscofile=paste0(dbxdir, '/qc/buscos.csv')
buscocsv=read_csv(buscofile) %>%
    select(-c(total))
busco=gather(buscocsv, key, value, -asm)

buscoplot=paste0(dbxdir, '/qc/buscos.pdf')
pdf(buscoplot, height=5, width=10)
plot=ggplot(busco, aes(x=asm, y=value, fill=key, colour=key, alpha=.5)) +
    geom_bar(width=.7, stat='identity') +
    coord_flip() +
    ggtitle('BUSCO') +
    xlab('Genome') +
    ylab('Number of buscos') +
    theme_bw()
print(plot)
dev.off()


##plot
countsdir=paste0(dbxdir, 'error/')
countfiles=list.files(countsdir, 'csv')
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
    dev.off()
    
    pdf(paste0(countsdir, name, '_table.pdf'))
    print(grid.table(counts[1:20,1:2], rows=NULL))
    dev.off()
}

mumdir='/kyber/Data/seqlab/sp_2019/nivar_r9/mummer_corr/allsnps/'
snpsfiles=list.files(mumdir, 'snps')
r9=foreach(i=snpsfiles, .combine=rbind) %dopar% {
    numcors=countLines(paste0(mumdir, i))
    return(data.frame(samp=i, cors=numcors))
}
r9$pore='r9'

mumdir='/kyber/Data/Nanopore/Analysis/190706_nivar_r10/mummer_corr/allsnps/'
snpsfiles=list.files(mumdir, 'snps')
r10=foreach(i=snpsfiles, .combine=rbind) %dopar% {
    numcors=countLines(paste0(mumdir, i))
    return(data.frame(samp=i, cors=numcors))
}
r10$pore='r10'

allcorr=rbind(r9, r10)

totalplot=paste0(dbxdir, '/qc/allerrors.pdf')
pdf(totalplot, height=8.5, width=11)
ggplot(allcorr, aes(x=pore, y=cors)) +
    geom_bar(aes(fill=samp, colour=samp, alpha=.8),width=.5, stat='identity', position='dodge') +
    ggtitle('Corrections Made') +
    theme_bw()
dev.off()


##plot busco
buscofile=paste0(dbxdir, '/qc/trans_buscos.csv')
buscocsv=read_csv(buscofile) %>%
    select(-c(total))
busco=gather(buscocsv, key, value, -asm)

buscoplot=paste0(dbxdir, '/qc/trans_buscos.pdf')
pdf(buscoplot, height=5, width=10)
plot=ggplot(busco, aes(x=asm, y=value, fill=key, colour=key, alpha=.5)) +
    geom_bar(width=.7, stat='identity') +
    coord_flip() +
    ggtitle('BUSCO') +
    xlab('Genome') +
    ylab('Number of buscos') +
    theme_bw()
print(plot)
dev.off()
