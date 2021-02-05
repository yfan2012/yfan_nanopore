library(tidyverse)

datadir='/uru/Data/Nanopore/projects/mdr'
krakendir='/uru/Data/Nanopore/projects/mdr/MDRstool_16/kraken'
rebasefile=file.path(datadir, 'refs', 'rebase_report.csv')

rebase=read_csv(rebasefile) %>%
    filter(!is.na(Specificity)) %>%
    filter(Gene=='M') %>%
    
motifcounts=rebase %>%
    group_by(Specificity) %>%
    summarise(count=n(), type=Type[1]) %>%
    arrange(-count)


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



test=rebase %>%
    rowwise() %>%
    filter(grepl('Bacteroides', Organism, fixed=TRUE))

motifsbyname <- function(rebase, name) {
    orgdata=rebase %>%
        filter(Organism==name)
    motifs=unique(orgdata$Specificity)
    return(motifs)
}

