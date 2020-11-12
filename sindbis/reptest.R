library(ggplot2)
library(tidyverse)
library(R.utils)
library(Biostrings)
library(ShortRead)

##read in cov files, basic cov plot
datadir='/dilithium/Data/Nanopore/sindbis'
dbxdir='~/Dropbox/timplab_data/sindbis/replicates/cov'

##annot_regions is a manually made regions file from the gb file in the dbox
gff=paste0(datadir, '/annot_regions.tsv')
regions=read_tsv(gff, col_names=c('prot', 'start', 'end')) %>%
    mutate(ymax=0) %>%
    mutate(colors=c('#EF9A9A', '#CE93D8', '#9FA8DA', '#81D4FA', '#80CBC4', '#C5E1A5', '#FFF59D', '#FFCC80', '#FFAB91'))

testsamps=c('reptest_sinv', 'reptest_sinvab')

getcov <- function(samps) {
    allcov=tibble(
        chr=as.character(),
        start=as.numeric(),
        end=as.numeric(),
        pos=as.numeric(),
        cov=as.numeric(),
        sample=as.character(),
        normcov=as.numeric())
    cnames=c('chr', 'start', 'end', 'pos', 'cov')
    for (samp in samps) {
        covfile=file.path(datadir, samp, 'cov', paste0(samp, '.primary.cov')) 
        cover=read_tsv(covfile, col_names=cnames) %>%
            mutate(sample=samp) %>%
            mutate(normcov=cov/sum(cov))
        
         allcov=bind_rows(allcov, cover)
    }
    return(allcov)
}

allcov=getcov(testsamps)
barheight=-.1*max(allcov$normcov)
regions$ymin=barheight

pdffile=file.path(dbxdir, 'reptest.pdf')
pdf(pdffile, w=18, height=5)
plot=ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot)) +
    geom_rect(show.legend=FALSE, alpha=.5) +
    geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
    scale_fill_manual(values=regions$colors, labels=regions$prot) +
    geom_line(data=allcov, inherit.aes=F, aes(x=pos, y=normcov, colour=sample), size=1)+
    scale_color_brewer(palette = "Set2") +
    xlab('Sindbis Genome Position') +
    ylab('Depth') +
    ggtitle(paste0('Replicates Test Coverage (normalized to total sindbis bases)')) +
    theme_bw()
print(plot)
dev.off()
