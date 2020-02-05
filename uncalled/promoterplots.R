library(doParallel)
library(tidyverse)
cl=makeCluster(12)
registerDoParallel(cl, cores=12)

uncalledfile='/uru/Data/Nanopore/projects/uncalled/methcalls/promoters_methcalls.tsv'
controlfile='/uru/Data/Nanopore/projects/uncalled/methcalls/promoters_ontwgs.tsv'
plotdir='~/Dropbox/timplab_data/uncalled/methplots/'

tsvcols=c('chr', 'methstart', 'methend', 'num_motifs', 'calledsites', 'calledmeth', 'methfreq', 'seq', 'prom_chr', 'meta', 'type', 'start', 'end', 'one', 'two', 'three', 'info')
uncalled=read_tsv(uncalledfile, col_names=tsvcols) %>%
    group_by(prom_chr, start, end) %>%
    summarise(meanmeth=mean(methfreq), numsites=length(methfreq)) %>%
    filter(numsites>20) %>%
    mutate(id=paste0(prom_chr,':', as.character(start),'-', as.character(end)))
control=read_tsv(controlfile, col_names=tsvcols) %>%
    group_by(prom_chr, start, end) %>%
    summarise(meanmeth=mean(methfreq), numsites=length(methfreq)) %>%
    filter(numsites>20) %>%
    mutate(id=paste0(prom_chr,':', as.character(start),'-', as.character(end)))

meth=uncalled %>%
    full_join(control, by='id') %>%
    drop_na()


pdf(paste0(plotdir, 'promoter_meth.pdf'))
p=ggplot(meth, aes(meanmeth.x, meanmeth.y)) +
    geom_bin2d() +
    ggtitle('Promoter Methylation') +
    xlab('UNCALLED') +
    ylab('ONT WGS') +
    theme_bw()
plot(p)
dev.off()



correlation=cor.test(meth$meanmeth.y, meth$meanmeth.x, method='pearson')
linmod=lm(meanmeth.x ~ meanmeth.y, data=meth)

##why is correlation only .91?
meth=meth %>%
    mutate(diff=meanmeth.x-meanmeth.y) %>%
    select(-c('prom_chr.y', 'start.y', 'end.y')) %>%
    arrange(-abs(diff))





uncalledfile='/uru/Data/Nanopore/projects/uncalled/methcalls/methcalls.tsv'
controlfile='/uru/Data/Nanopore/projects/uncalled/methcalls/control_methcalls.tsv'

tsvcols=c('chr', 'strand', 'start', 'end', 'readname', 'llr', 'llm', 'llu', 'numstrands', 'nummotifs', 'seq')

uncalledraw=read_tsv(uncalledfile, col_names=tsvcols, skip=1) %>%
    select('llr') %>%
    mutate(samp='uncalled')
controlraw=read_tsv(controlfile, col_names=tsvcols, skip=1) %>%
    select('llr') %>%
    mutate(samp='control')

ratios=rbind(controlraw, uncalledraw)
pdf(paste0(plotdir, 'nanopolish_llr.pdf'))
p=ggplot(ratios, aes(x=llr, colour=samp, fill=samp, alpha=.2)) +
    geom_density() +
    xlim(-15, 15) +
    ggtitle('Nanopolish Loglik Ratio Distributions') +
    theme_bw()
plot(p)
dev.off()

