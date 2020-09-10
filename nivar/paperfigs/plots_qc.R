library(tidyverse)
library(tidyr)
library(gridExtra)
library(Biostrings)
library(R.utils)

dbxdir='~/Dropbox/yfan/nivar/'
datadir='/uru/Data/Nanopore/projects/nivar/'
asmfile=paste0(datadir,'pilon/r9_pilon/nivar_r9.pilon_bwa.6.fasta')

##plot busco
cnames=c('buscoid', 'status', 'contig', 'start', 'end', 'score', 'length')
asmbuscofile=paste0(datadir, 'busco/run_r9_pilon/full_table_r9_pilon.tsv')
asmbusco=read_tsv(asmbuscofile, comment='#', col_names=cnames) %>%
    group_by(status) %>%
    summarise(sum=n()) %>%
    mutate(asm='asm')
refbuscofile=paste0(datadir, 'busco/run_ref_busco/full_table_ref_busco.tsv')
refbusco=read_tsv(refbuscofile, comment='#', col_names=cnames) %>%
    group_by(status) %>%
    summarise(sum=n()) %>%
    mutate(asm='ref')
busco=rbind(refbusco, asmbusco)

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


##check telomeres
asm=readDNAStringSet(asmfile)
