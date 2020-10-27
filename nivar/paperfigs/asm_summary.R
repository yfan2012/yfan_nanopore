library(tidyverse)
library(argparse)
library(Biostrings)

##misc paper figs
dbxdir='~/Dropbox/yfan/nivar/paperfigs/raw/'
datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/'

reffile='/uru/Data/Nanopore/projects/nivar/reference/candida_nivariensis.fa'
asmfile=paste0(datadir,'assembly/nivar.contigs.fasta')
finfile=paste0(datadir, '/assembly_final/nivar.final.fasta')

##asm contig lengths compared to ref length
fin=readDNAStringSet(finfile)
fininfo=tibble(lens=sort(width(fin), decreasing=TRUE)) %>%
    mutate(rank=row_number()) %>%
    mutate(cumulative=cumsum(lens)) %>%
    mutate(perc=cumulative/max(cumulative)) %>%
    mutate(samp='fin')

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

##busco plots and comparisons
cnames=c('buscoid', 'status', 'contig', 'start', 'end', 'score', 'length')
scabuscofile=paste0(datadir, 'busco/run_nivar_scaffold/full_table_nivar_scaffold.tsv')
scabusco=read_tsv(scabuscofile, comment='#', col_names=cnames) %>%
    group_by(status) %>%
    summarise(sum=n()) %>%
    mutate(asm='asm')
refbuscofile=paste0(datadir, 'busco/run_ref/full_table_ref.tsv')
refbusco=read_tsv(refbuscofile, comment='#', col_names=cnames) %>%
    group_by(status) %>%
    summarise(sum=n()) %>%
    mutate(asm='ref')
busco=rbind(refbusco, scabusco)

buscoplot=paste0(dbxdir, 'paperfigs/raw/buscos.pdf')
pdf(buscoplot, height=4, width=13)
plot=ggplot(busco, aes(x=asm, y=sum, fill=status, colour=status, alpha=.5)) +
    geom_bar(width=.7, stat='identity') +
    coord_flip() +
    ggtitle('BUSCO') +
    xlab('Genome') +
    ylab('Number of buscos') +
    theme_bw()
print(plot)
dev.off()

if (FALSE) {
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
}


