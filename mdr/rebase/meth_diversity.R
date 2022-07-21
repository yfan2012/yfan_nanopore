library(tidyverse)
library(RColorBrewer)

srcdir='~/Code/yfan_nanopore/mdr/rebase'
datadir='/uru/Data/Nanopore/projects/mdr'
krakendir='/uru/Data/Nanopore/projects/mdr/MDRstool_16/kraken'
dbxdir='~/gdrive/mdr/rebase'
rebasefile=file.path(datadir, 'refs', 'rebase_report.csv')

motifsbyname <- function(rebase, name) {
    orgdata=rebase %>%
        filter(Organism==name)
    motifs=unique(orgdata$Specificity)
    return(motifs)
}

rebase=read_csv(rebasefile) %>%
    filter(!is.na(Specificity)) %>%
    filter(Gene=='M')
rebasefiltcsv=file.path(datadir, 'refs', 'rebase_filtered.csv')
write_csv(rebase, rebasefiltcsv)

common=sort(table(rebase$Specificity), decreasing=TRUE)

freq=sort(table(rebase$Specificity), decreasing=T)


motifcounts=rebase %>%
    group_by(Specificity) %>%
    summarise(count=n(), type=Type[1]) %>%
    arrange(-count) %>%
    filter(type!='I') %>%
    mutate(perc=count/sum(count))
top10file=file.path(srcdir, 'barcodes10.txt')
write_tsv(motifcounts[1:10,1], top10file, col_names=FALSE)
top15file=file.path(srcdir, 'barcodes15.txt')
write_tsv(motifcounts[1:15,1], top15file, col_names=FALSE)
top20file=file.path(srcdir, 'barcodes20.txt')
write_tsv(motifcounts[1:21,1], top20file, col_names=FALSE)
top50file=file.path(srcdir, 'barcodes50.txt')
write_tsv(motifcounts[1:50,1], top50file, col_names=FALSE)

countsfile=file.path(dbxdir, 'motif_abundance.pdf')
pdf(countsfile, h=9, w=15)
motiforder=motifcounts$Specificity[1:20]
ggplot(motifcounts[1:20,], aes(x=Specificity, y=count, colour=type, fill=type, alpha=.5)) +
    geom_bar(stat='identity') +
    ggtitle('Motif abundances') +
    xlab('Motif') +
    ylab('Abundance') +
    scale_x_discrete(limits = motiforder) +
    scale_fill_brewer(palette = 'Set2') +
    scale_colour_brewer(palette = 'Set2') +
    theme_bw()
ggplot(motifcounts[1:20,], aes(x=Specificity, y=perc, colour=type, fill=type, alpha=.5)) +
    geom_bar(stat='identity') +
    ggtitle('Motif abundances') +
    xlab('Motif') +
    ylab('Abundance') +
    scale_x_discrete(limits = motiforder) +
    scale_fill_brewer(palette = 'Set2') +
    scale_colour_brewer(palette = 'Set2') +
    theme_bw()
dev.off()


datasets=c('native', 'pcr', 'phase', 'shotgun')
krakeninfo=tibble(percent=as.numeric(),
                  numcov=as.integer(),
                  numclass=as.integer(),
                  rank=as.character(),
                  id=as.integer(),
                  name=as.character())
reportcols=c('percent', 'numcov', 'numclass', 'rank', 'id', 'name')
for (i in datasets) {
    reportfile=file.path(krakendir, paste0(i, '.report.top40.txt'))
    setreport=read_tsv(reportfile, col_names=reportcols)
    krakeninfo=bind_rows(krakeninfo, setreport)
}

speciesinfo=krakeninfo %>%
    group_by(name) %>%
    summarise(name=name[1], rank=rank[1], prevalance=sum(percent)) %>%
    arrange(-prevalance) %>%
    filter(name!='Homo sapiens')


species_motifs=rebase %>%
    rowwise() %>%
    filter(Organism %in% speciesinfo$name) %>%
    filter(Type=='II', Gene=='M')
motifs=unique(species_motifs$Specificity)


genomelist=tibble(acc=unique(rebase$GenBank))
methdbdir='/atium/Data/ref/mtases'
listfile=file.path(methdbdir, 'genomelist.txt')
write_tsv(genomelist, listfile, col_names=FALSE)
