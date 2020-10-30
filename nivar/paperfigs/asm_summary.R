library(tidyverse)
library(RColorBrewer)
library(ggforce)
library(argparse)
library(Biostrings)
library(RIdeogram)
source('~/Code/yfan_nanopore/nivar/paperfigs/genome.R')

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
asms=rbind(fininfo, refinfo)

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





##genome busco plots
cnames=c('buscoid', 'status', 'contig', 'start', 'end', 'score', 'length', 'url', 'description')
buscos=c('asm', 'ref', 'gla', 'cer', 'alb')
allsum=tibble(
    status=as.character(),
    sum=as.numeric(),
    asm=as.character())
for (i in buscos) {
    buscofile=file.path(datadir, 'busco', i, 'run_saccharomycetes_odb10/full_table.tsv')
    sum=read_tsv(buscofile, comment='#', col_names=cnames) %>%
        group_by(status) %>%
        summarise(sum=n()) %>%
        mutate(asm=i)
    allsum=rbind(allsum, sum)
}
allsum=allsum %>%
    mutate(alpha=case_when(asm=='asm' | asm=='ref' ~ 1, T ~ 0)) %>%
    mutate(order=case_when(asm=='asm' ~ 5,
                           asm=='ref' ~ 4,
                           asm=='gla' ~ 3,
                           asm=='cer' ~ 2,
                           asm=='alb' ~ 1))
buscoplot=paste0(dbxdir, 'buscos.pdf')
pdf(buscoplot, height=8, width=19)
plot=ggplot(allsum, aes(x=status, y=sum, fill=asm, group=order, alpha=alpha)) +
    geom_bar(position='dodge', stat='identity', width=.7) +
    scale_alpha_continuous(range=c(.25,.8)) +
    scale_fill_brewer(palette = "Set2") +
    ggtitle('BUSCO') +
    xlab('BUSCO Status') +
    ylab('Number of BUSCOs') +
    theme_bw() +
    facet_zoom(ylim = c(0, 100))
print(plot)
dev.off()

##check if buscos are in common
refbuscofile=file.path(datadir, 'busco', 'ref', 'run_saccharomycetes_odb10/full_table.tsv')
refmiss=read_tsv(refbuscofile, comment='#', col_names=cnames) %>%
    filter(status=='Missing')
finbuscofile=file.path(datadir, 'busco', 'asm', 'run_saccharomycetes_odb10/full_table.tsv')
finmiss=read_tsv(finbuscofile, comment='#', col_names=cnames) %>%
    filter(status=='Missing')
    

##transcriptome busco
transbuscos=c('asm', 'gla', 'cer', 'alb')
alltrans=tibble(
    status=as.character(),
    sum=as.numeric(),
    asm=as.character())
for (i in transbuscos) {
    buscofile=file.path(datadir, 'busco', paste0('trans', i), 'run_saccharomycetes_odb10/full_table.tsv')
    sum=read_tsv(buscofile, comment='#', col_names=cnames) %>%
        group_by(status) %>%
        summarise(sum=n()) %>%
        mutate(asm=i)
    alltrans=rbind(alltrans, sum)
}
alltrans=alltrans %>%
    mutate(alpha=case_when(asm=='asm' ~ 1, T ~ 0)) %>%
    mutate(order=case_when(asm=='asm' ~ 5,
                           asm=='gla' ~ 3,
                           asm=='cer' ~ 2,
                           asm=='alb' ~ 1))
buscoplot_trans=paste0(dbxdir, 'buscos_transcriptome.pdf')
pdf(buscoplot_trans, height=8, width=19)
plot=ggplot(alltrans, aes(x=status, y=sum, fill=asm, group=order, alpha=alpha)) +
    geom_bar(position='dodge', stat='identity', width=.6) +
    scale_alpha_continuous(range=c(.25,.8)) +
    scale_fill_brewer(palette = "Set2") +
    ggtitle('BUSCO') +
    xlab('BUSCO Status') +
    ylab('Number of BUSCOs') +
    facet_zoom(ylim = c(0, 100)) +
    theme_bw()
print(plot)
dev.off()


##ideogram
ideodata=tibble(
    Chr=names(fin)[1:15],
    Start=rep(0, times=15),
    End=width(fin)[1:15])
ideowindows=tibble(
    Chr=as.character(),
    Start=as.numeric(),
    End=as.numeric())
for (i in 1:dim(ideodata)[1]) {
    tigranges=tibble(
        Start=seq(1, ideodata$End[i], 50000), 
        End=c(seq(20001, ideodata$End[i], 20000), ideodata$End[i])) %>%
        mutate(Chr=ideodata$Chr[i])
    ideowindows=rbind(ideowindows, tigranges)
}
ideowindows=ideowindows %>%
    mutate(numfwd=str_count(substr(as.character(fin[Chr]), Start, End), fwdtelo)) %>%
    mutate(numrev=str_count(substr(as.character(fin[Chr]), Start, End), revtelo)) %>%
    mutate(numtelo=numfwd+numrev)
ideotelo=tibble(
    Chr=ideowindows$Chr,
    Start=ideowindows$Start,
    End=ideowindows$End,
    Value=ideowindows$numtelo) %>%
    mutate(color=case_when(Value > '1' ~ 'ff0000', T ~ '000000'))

ideosmall=tibble(
    Chr=as.character(),
    Start=as.numeric(),
    End=as.numeric())
for (i in 1:dim(ideodata)[1]) {
    tigranges=tibble(
        Start=seq(1, ideodata$End[i], 500), 
        End=c(seq(501, ideodata$End[i], 500), ideodata$End[i])) %>%
        mutate(Chr=ideodata$Chr[i])
    ideosmall=rbind(ideosmall, tigranges)
}
ideosmall=ideosmall %>%
    mutate(numA=str_count(substr(as.character(fin[Chr]), Start, End), 'A')) %>%
    mutate(numT=str_count(substr(as.character(fin[Chr]), Start, End), 'T')) %>%
    mutate(total=numA+numT)
overlay=ideosmall %>%
    select(Chr, Start, End) %>%
    mutate(Value=1/(1+exp(-(ideosmall$total/max(ideosmall$total)))))

drawideo=ideogram(karyotype=as.data.frame(ideodata),
                  overlaid=as.data.frame(overlay), 
                  label=(as.data.frame(ideotelo)),
                  label_type='polygon',
                  colorset1=c("#ff0000", "#ffffff", "#0000ff"),
                  output=file.path(dbxdir,'ideogram.svg'))
                 


telocheck <- function(asmfile, fwdtelo, revtelo) {
    asm=readDNAStringSet(asmfile)
    seqs=as.character(asm)

    teloinfo=tibble(chr=names(asm)) %>%
        mutate(length=width(asm[chr])) %>%
        mutate(fwd=str_count(as.character(asm[chr]), fwdtelo)) %>%
        mutate(rev=str_count(as.character(asm[chr]), revtelo))
}
