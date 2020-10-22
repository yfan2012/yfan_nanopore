library(tidyverse)
library(argparse)
library(Biostrings)

##misc paper figs
dbxdir='~/Dropbox/yfan/nivar/paperfigs/raw/'
datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/'
reffile='/uru/Data/Nanopore/projects/nivar/reference/candida_nivariensis.fa'

##asm contig lengths compared to ref length
asmfile=paste0(datadir, '/assembly_final/nivar.final.fasta')
asm=readDNAStringSet(asmfile)
asminfo=tibble(lens=sort(width(asm), decreasing=TRUE)) %>%
    mutate(rank=row_number()) %>%
    mutate(cumulative=cumsum(lens)) %>%
    mutate(perc=cumulative/max(cumulative)) %>%
    mutate(samp='asm')
ref=readDNAStringSet(reffile)
refinfo=tibble(lens=sort(width(ref), decreasing=TRUE)) %>%
    mutate(rank=row_number()) %>%
    mutate(cumulative=cumsum(lens)) %>%
    mutate(perc=cumulative/max(cumulative)) %>%
    mutate(samp='ref')
asms=rbind(asminfo, refinfo)

pdffile=paste0(dbxdir, 'contig_lengths.pdf')
pdf(pdffile, w=12, h=7)
ggplot(asms, aes(x=rank, y=perc, colour=samp)) +
    geom_step(size=1) +
    ggtitle('Cumulative Percentage of Assembly Length')  +
    xlim(0,50) +
    xlab('Contig Rank') +
    ylab('Percent of Assembly Length') +
    theme_bw()
ggplot(asms, aes(x=rank, y=perc, colour=samp)) +
    geom_step(size=1) +
    ggtitle('Cumulative Percentage of Assembly Length')  +
    xlab('Contig Rank') +
    ylab('Percent of Assembly Length') +
    theme_bw()
dev.off()

    
