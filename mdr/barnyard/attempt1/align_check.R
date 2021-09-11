library(tidyverse)

datadir='/mithril/Data/Nanopore/projects/methbin/barnyard/align'
dbxdir='~/gdrive/mdr/barnyard'

pafcols=c('qname', 'qlen', 'qstart', 'qend', 'strand', 'rname', 'rlen', 'rstart', 'rend', 'num_match', 'alen', 'mapq', 'tp', 'cm', 's1', 's2', 'dv', 'rl')

aligninfo=tibble()
for (i in c(1,2,3,4)) {
    prefix=paste0('210730_mdr_barnyard_mix', as.character(i))
    paffile=file.path(datadir, paste0(prefix, '.paf'))

    pafinfo=read_tsv(paffile, col_names=pafcols) %>%
        filter(mapq>40) %>%
        select(qlen, alen, rname) %>%
        mutate(samp=paste0('mix',as.character(i)))
    aligninfo=bind_rows(aligninfo, pafinfo)
}

##plot lengths
lengthsfile=file.path(dbxdir, '210730_runs_lengths.pdf')
pdf(lengthsfile, w=13, h=7)
plot=ggplot(aligninfo, aes(x=alen, colour=rname, fill=rname, alpha=.2)) +
    geom_density() +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    facet_wrap(~samp) +
    theme_bw()
print(plot)
dev.off()



                      
