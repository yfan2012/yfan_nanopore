library(doParallel)
library(tidyverse)
cl=makeCluster(12)
registerDoParallel(cl, cores=12)

datadir='/uru/Data/Nanopore/projects/uncalled/'

##check uncalled vs control
uncalledfile=paste0(datadir,'methcalls/promoters_methcalls.tsv')
controlfile=paste0(datadir, 'methcalls/promoters_ontwgs.tsv')
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

meth=meth %>%
    mutate(diff=meanmeth.x-meanmeth.y) %>%
    select(-c('prom_chr.y', 'start.y', 'end.y')) %>%
    arrange(-abs(diff))





##check uncalled vs encode bisulfite
wgbs1file=paste0(datadir,'bisulfite/promoters_ENCFF279HCL.tsv')
wgbs2file=paste0(datadir,'bisulfite/promoters_ENCFF835NTC.tsv')

cols=c('chr', 'methstart', 'methend', 'item', 'score', 'strand','startcod', 'endcod', 'color', 'cov', 'methfreq', 'prom_chr', 'meta', 'type', 'start', 'end', 'one', 'two', 'three', 'info')
wgbs1=read_tsv(wgbs1file, col_names=cols) %>%
    group_by(prom_chr, start, end) %>%
    summarise(meanmeth=mean(methfreq), meancov=mean(cov), numsites=length(methfreq)) %>%
    filter(numsites>20) %>%
    mutate(id=paste0(prom_chr,':', as.character(start),'-', as.character(end)))
wgbs2=read_tsv(wgbs2file, col_names=cols) %>%
    group_by(prom_chr, start, end) %>%
    summarise(meanmeth=mean(methfreq), meancov=mean(cov), numsites=length(methfreq)) %>%
    filter(numsites>20) %>%
    mutate(id=paste0(prom_chr,':', as.character(start),'-', as.character(end)))

wgbsnp1=uncalled %>%
    full_join(wgbs1, by='id') %>%
    drop_na()
wgbsnp2=uncalled %>%
    full_join(wgbs2, by='id') %>%
    drop_na()
wgbsonly=wgbs1 %>%
    full_join(wgbs2, by='id') %>%
    drop_na() %>%
    filter(id %in% uncalled$id)


pdf(paste0(plotdir, 'promoter_meth_bisulfite.pdf'))
p=ggplot(wgbsnp1, aes(meanmeth.x, meanmeth.y)) +
    geom_bin2d() +
    ggtitle('Promoter Methylation') +
    xlab('UNCALLED') +
    ylab('WGBS1') +
    theme_bw()
plot(p)
q=ggplot(wgbsnp2, aes(meanmeth.x, meanmeth.y)) +
    geom_bin2d() +
    ggtitle('Promoter Methylation') +
    xlab('UNCALLED') +
    ylab('WGBS2') +
    theme_bw()
plot(q)
r=ggplot(wgbsonly, aes(meanmeth.x, meanmeth.y)) +
    geom_bin2d() +
    ggtitle('Promoter Methylation') +
    xlab('WGBS1') +
    ylab('WGBS2') +
    theme_bw()
plot(r)
dev.off()

correlation=cor.test(wgbsnp1$meanmeth.y, wgbsnp1$meanmeth.x, method='pearson')
correlation=cor.test(wgbsnp2$meanmeth.y, wgbsnp2$meanmeth.x, method='pearson')
correlation=cor.test(wgbsonly$meanmeth.y, wgbsonly$meanmeth.x, method='pearson')



##check log lik ratio distributions
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

