library(tidyverse)
library(umap)
library(RColorBrewer)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/figs'


##chrom to bin info
chrombinsfile='/mithril/Data/Nanopore/projects/methbin/paperfigs/contig_level/tigs2bins.tsv'
chrombins=read_tsv(chrombinsfile)

##read plas stuff
plasinfofile=file.path('~/gdrive/mdr/paperfigs/contig_level/nearestplas_cov.csv')
plasinfo=read_csv(plasinfofile) %>%
    filter(mumbin!='unknown') %>%
    filter(nummoitfs>=10)


##set up frequency info
exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))
freqs=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec) %>%
    spread(key=motif, value=freq)

plasfreqs=freqs %>%
    filter(chrom %in% plasinfo$chroms) %>%
    filter(complete.cases(.))

nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])

methchroms=methfreq %>%
    filter(chrom %in% keepchroms)
chrominfo=methchroms %>%
    spread(key=motif, value=freq) %>%
    filter(chrom %in% chrombins$rname)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom


##umap
chromumap=umap(matchrominfo)
chromplot=tibble(chroms=rownames(chromumap$layout)) %>%
    mutate(x=chromumap$layout[,1]) %>%
    mutate(y=chromumap$layout[,2]) %>%
    rowwise() %>%
    mutate(bin=chrombins$bin[chrombins$rname==chroms]) %>%
    mutate(tiglen=chrombins$rlen[chrombins$rname==chroms]) %>%
    mutate(type=case_when(!chroms %in% plasinfo$chroms ~ 'chr',
                          TRUE ~ 'plas'))


##plot
mycolors=read_csv('colors.csv')
cols=as.vector(mycolors$color)
names(cols)=mycolors$bin

umappdf=file.path(dbxdir, 'umap.pdf')
pdf(umappdf, h=8, w=11)
plot=ggplot(chromplot %>% filter(type=='chr'), aes(x=x, y=y, colour=bin, size=tiglen, alpha=.1)) +
    geom_point() +
    scale_colour_manual(values=cols) +
    ggtitle('contig umap') +
    theme_bw()
mainplot=plot +
    geom_point(data=chromplot %>% filter(type=='plas'), aes(x=x, y=y, size=tiglen), shape=4, inherit.aes=FALSE)
print(mainplot)
dev.off()
    


