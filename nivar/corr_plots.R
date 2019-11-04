library(ggplot2)
library(tidyverse)
library(tidyr)
library(gridExtra)
library(doParallel)
library(R.utils)

dbxdir='~/Dropbox/yfan/nivar/'
countsdir=paste0(dbxdir, 'motif_enrich/counts/')

##make tables of kmer populations
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


##make venn diagrams comparing kmer distributions
corrs=c('pilon', 'racon', 'freebayes')

allranks=tibble(
    topmers = numeric(),
    common = numeric(),
    common_perc = numeric(),
    correction = character(), 
    r9perc = numeric(),
    r10perc = numeric()
)

corrpot=tibble(
    seq = character(),
    diffrank = numeric(),
    diffperc = numeric(),
    correction = character()
)
    
for (i in corrs) {
    r9=paste0(countsdir, '/nivar_r9_', i, '_raw.csv')
    r10=paste0(countsdir, '/nivar_r10_', i, '_raw.csv')

    r9counts=read_csv(r9, col_names=c('seq', 'freq')) %>%
        arrange(-freq) %>%
        mutate(perc_errors=freq/sum(freq)) %>%
        mutate(rank=1:length(freq)) %>%
        mutate(totalerr=1-cumsum(perc_errors))
    
    r10counts=read_csv(r10, col_names=c('seq', 'freq')) %>%
        arrange(-freq) %>%
        mutate(perc_errors=freq/sum(freq)) %>%
        mutate(rank=1:length(freq)) %>%
        mutate(totalerr=1-cumsum(perc_errors))

    pot9=r9counts %>%
        arrange(seq)
    
    pot10=r10counts %>%
        arrange(seq) %>%
        mutate(diffrank=abs(rank-pot9$rank)) %>%
        mutate(diffperc=abs(perc_errors-pot9$perc_errors)) %>%
        mutate(correction=i)

    corrpot=bind_rows(pot10[,c(1,6,7,8)])
    
    ranks=tibble(topmers=1:4096)
    ranks=ranks %>%
        rowwise() %>%
        mutate(common=sum(r9counts$seq[1:topmers] %in% r10counts$seq[1:topmers])) %>%
        mutate(common_perc=common/topmers) %>%
        mutate(correction=i) %>%
        add_column(r9perc=r9counts$totalerr[1:4096]) %>%
        add_column(r10perc=r10counts$totalerr[1:4096])
    
    allranks=bind_rows(allranks, ranks)
}    

pdf(paste0(dbxdir, 'motif_enrich/common_kmers.pdf'), height=8.5, width=11)
plot=ggplot(allranks, aes(x=topmers, y=common_perc, colour=correction)) +
    geom_line() +
    geom_line(linetype='dotted', aes(x=topmers, y=r9perc)) +
    geom_line(linetype='twodash',aes(x=topmers, y=r10perc)) +
    ggtitle('Erroneous 6-mers  in common between R9 and R10 assemblies') +
    xlab('6-mer rank') +
    ylab('Percent 6-mers in common') +
    ylim(0,1) +
    theme_bw()
print(plot)
dev.off()

pdf(paste0(dbxdir, 'motif_enrich/correction_potential.pdf'), height=8.5, width=11)
rankplot=ggplot(corrpot, aes(x=diffrank)) +
    geom_histogram() +
    ggtitle('Differences in 6-mer rank wrt error') +
    xlab('Difference in 6mer rank between r9 and r10') +
    ylab('Count') +
    theme_bw()
print(rankplot)
percplot=ggplot(corrpot, aes(x=diffperc)) +
    geom_histogram() +
    ggtitle('Differences in how much error each 6-mer accounts for ') +
    xlab('Difference in percent error accounted for in r9 vs r10') +
    ylab('Count') +
    theme_bw()
print(percplot)
dev.off()
