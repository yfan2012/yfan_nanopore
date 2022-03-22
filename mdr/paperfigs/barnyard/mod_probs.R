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


modprobdistpdf=file.path(dbxdir, 'methprobs.pdf')
pdf(modprobdistpdf, w=11, h=8)
damplot=ggplot(daminfo, aes(x=modprob, colour=samp, fill=samp, alpha=.1)) +
    geom_density() +
    ggtitle('dam') +
    facet_wrap(~pos) +
    theme_bw()
print(damplot)
dcmplot=ggplot(dcminfo, aes(x=modprob, colour=samp, fill=samp, alpha=.1)) +
    geom_density() +
    ggtitle('dcm') +
    facet_wrap(~pos) +
    theme_bw()
print(dcmplot)
dev.off()
