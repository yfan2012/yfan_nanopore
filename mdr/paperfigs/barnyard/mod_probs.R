library(tidyverse)

datadir='/mithril/Data/Nanopore/projects/methbin/barnyard/strains/modcheck'
dbxdir='~/gdrive/mdr/barnyard'
files=list.files(datadir)


methprobs=NULL
modcols=c('readname', 'chrom', 'strand', 'pos', 'logmodprob', 'logumodprob', 'modbase')
for (i in files) {
    chr=strsplit(i, '_')[[1]][1]
    motif=strsplit(strsplit(i, '_')[[1]][2], '[.]')[[1]][1]
    modfile=file.path(datadir, i)

    modinfo=read_tsv(modfile, col_names=modcols) %>%
        mutate(modprob=exp(logmodprob)) %>%
        mutate(motif=motif) %>%
        mutate(samp=chr)

    methprobs=bind_rows(methprobs, modinfo)
}

daminfo=methprobs %>%
    filter(motif=='dam')
dcminfo=methprobs %>%
    filter(motif=='dcm')

dbxdir='~/gdrive/mdr/paperfigs/figs'

modprobdistpdf=file.path(dbxdir, 'barnyard_methprobs.pdf')
pdf(modprobdistpdf, w=11, h=8)
damplot=ggplot(daminfo, aes(x=modprob, colour=samp, fill=samp, alpha=.1)) +
    geom_density() +
    ggtitle('dam') +    
    facet_wrap(~pos) +
    geom_vline(xintercept=.2) +
    geom_vline(xintercept=.8) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
print(damplot)
dcmplot=ggplot(dcminfo, aes(x=modprob, colour=samp, fill=samp, alpha=.1)) +
    geom_density() +
    ggtitle('dcm') +
    facet_wrap(~pos) +
    geom_vline(xintercept=.2) +
    geom_vline(xintercept=.8) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
print(dcmplot)
dev.off()


dam.midcounts=daminfo %>%
    group_by(chrom, pos, samp) %>%
    summarise(middle=sum(modprob<.8 & modprob>.2)/n()) %>%
    mutate(motif='GATC')
dcm.midcounts=dcminfo %>%
    group_by(chrom, pos, samp) %>%
    summarise(middle=sum(modprob<.8 & modprob>.2)/n()) %>%
    mutate(motif='CCWGG')

midcounts=rbind(dam.midcounts, dcm.midcounts)
midcountscsv=file.path(dbxdir, 'barnyard_methprobs.csv')
write_csv(midcounts, midcountscsv)

    
