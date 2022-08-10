library(tidyverse)

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'


####assign contigs to bins based on the mummer
mumfile=file.path(datadir, 'clin_mummer', 'asm_hiC.1coords')
mumcols=c('rstart', 'rend', 'qstart', 'qend', 'ralen', 'qalen','ident', 'rlen', 'qlen', 'rcov', 'qcov', 'rname', 'qname')
mum=read_tsv(mumfile, col_names=mumcols) %>%
    rowwise() %>%
    mutate(bin=strsplit(qname, '.', fixed=TRUE)[[1]][1]) %>%
    group_by(rname, bin, rlen) %>%
    summarise(total_rcov=sum(rcov)) %>%
    filter(total_rcov>20)

mumkey=mum %>%
    group_by(rname) %>%
    filter(n()==1)

mummultiple=mum %>%
    group_by(rname, rlen) %>%
    filter(n()>1) %>%
    summarise(all_bins=paste0(bin, collapse=','), all_rcov=paste0(total_rcov, collapse=','))

mumkeyfile=file.path(datadir, 'tigs2bins.tsv')
write_tsv(mumkey, mumkeyfile)
multikeyfile=file.path(datadir, 'tigs2bins_multi.tsv')
write_tsv(mummultiple, multikeyfile)

##include mcoords file stuff
allOverlaps <- function(x) {
    ##x is a block of mummer output grouped by rname and bin                                                                                                                                                                                                                  
    ##returns one row tibble with [rname, bin, total cov]                                                                                                                                                                                                                     
    coveredpos=c()
    for (i in 1:dim(x)[1]) {
        coveredpos=c(coveredpos, seq(x$rstart[i], x$rend[i], by=1))
    }
    totcov=length(unique(coveredpos))
    return(tibble(rname=x$rname[1], bin=x$bin[1], rlen=x$rlen[1], total_rcov=100*totcov/x$rlen[1]))
    ##return(totcov/x$rlen[1])                                                                                                                                                                                                                                                
}

mumfilem=file.path(datadir, 'clin_mummer', 'asm_hiC.mcoords')
inclusive=read_tsv(mumfilem, col_names=mumcols) %>%
    rowwise() %>%
    mutate(bin=strsplit(qname, '.', fixed=TRUE)[[1]][1]) %>%
    ungroup() %>%
    group_by(rname, bin) %>%
    do(allOverlaps(.)) %>%
    filter(total_rcov>20) %>%
    filter(!rname %in% mum$rname)

mum.inclusive=rbind(mum, inclusive)
mumkey.inclusive=mum.inclusive %>%
    group_by(rname) %>%
    filter(n()==1)

mummultiple.inclusive=mum.inclusive %>%
    group_by(rname, rlen) %>%
    filter(n()>1) %>%
    summarise(all_bins=paste0(bin, collapse=','), all_rcov=paste0(total_rcov, collapse=','))

mumkeyfile=file.path(datadir, 'tigs2bins_inclusive.tsv')
write_tsv(mumkey.inclusive, mumkeyfile)
multikeyfile=file.path(datadir, 'tigs2bins_multi_inclusive.tsv')
write_tsv(mummultiple.inclusive, multikeyfile)





CATfile=file.path(datadir, 'contig_id/CAT/200708_mdr_stool16native_nohuman.CAT.names_official.txt')
phyloranks=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
CATcols=c('tig', 'classification', 'reason', 'lineage', 'lineage_scores', phyloranks)
CAT=read_tsv(CATfile)
names(CAT)=CATcols

BATfile=file.path(projdir, 'mdr/hiC/bin_id/BAT_single/200708_mdr_stool16native.BAT.names_official.txt')
BAT=read_tsv(BATfile)
names(BAT)=CATcols


tiginfo=NULL
for (i in 1:dim(mum)[1]) {
    tig=mum$rname[i]
    bin=mum$bin[i]

    infobin=BAT[BAT$tig==paste0(bin, '.fasta'),][-(1:5)] %>% slice(1) %>% unlist(., use.names=FALSE)
    infotig=CAT[CAT$tig==tig,][-(1:5)] %>% slice(1) %>% unlist(., use.names=FALSE)
    allinfo=tibble(bin=infobin, tig=infotig) %>%
        rowwise() %>%
        mutate(bin=strsplit(bin, ':', fixed=TRUE)[[1]][1]) %>%
        mutate(tig=strsplit(tig, ':', fixed=TRUE)[[1]][1])
    allinfo[is.na(allinfo)]='not applicable'
    
    binindex=sum(!allinfo$bin=='no support')
    tigindex=sum(!allinfo$tig=='no support')
    
    binleaf=allinfo$bin[binindex]
    tigleaf=allinfo$tig[tigindex]

    
    if (length(tigleaf)==0){
        tigleaf='none'
    }
    if (length(binleaf)==0){
        binleaf='none'
    }
    
    tiginfo=bind_rows(tiginfo, tibble(tig=tig, bin=bin, tigleaf=tigleaf, binleaf=binleaf))
}
names(mum)=c('tig', 'bin', 'rlen', 'total_rcov')

fullinfo=full_join(mum, tiginfo, by=c('tig', 'bin')) %>%
    arrange(bin)

speciesfile=file.path(datadir, 'tigbin_species.tsv')
write_tsv(fullinfo, speciesfile)




tiginfo.inclusive=NULL
for (i in 1:dim(mum.inclusive)[1]) {
    tig=mum.inclusive$rname[i]
    bin=mum.inclusive$bin[i]

    infobin=BAT[BAT$tig==paste0(bin, '.fasta'),][-(1:5)] %>% slice(1) %>% unlist(., use.names=FALSE)
    infotig=CAT[CAT$tig==tig,][-(1:5)] %>% slice(1) %>% unlist(., use.names=FALSE)
    allinfo=tibble(bin=infobin, tig=infotig) %>%
        rowwise() %>%
        mutate(bin=strsplit(bin, ':', fixed=TRUE)[[1]][1]) %>%
        mutate(tig=strsplit(tig, ':', fixed=TRUE)[[1]][1])
    allinfo[is.na(allinfo)]='not applicable'

    binindex=sum(!allinfo$bin=='no support')
    tigindex=sum(!allinfo$tig=='no support')

    binleaf=allinfo$bin[binindex]
    tigleaf=allinfo$tig[tigindex]


    if (length(tigleaf)==0){
        tigleaf='none'
    }
    if (length(binleaf)==0){
        binleaf='none'
    }

    tiginfo.inclusive=bind_rows(tiginfo.inclusive, tibble(tig=tig, bin=bin, tigleaf=tigleaf, binleaf=binleaf))
}

names(mum.inclusive)=c('tig', 'bin', 'rlen', 'total_rcov')
fullinfo.inclusive=full_join(mum.inclusive, tiginfo.inclusive, by=c('tig', 'bin')) %>%
    arrange(bin)

speciesfile.inclusive=file.path(datadir, 'tigbin_species_inclusive.tsv')
write_tsv(fullinfo.inclusive, speciesfile.inclusive)



##try with kraken
krakenfile=file.path(datadir, 'contig_id/kraken_standard/200708_mdr_stool16native_nohuman.out.txt')
kraken_cols=c('class', 'tig', 'tigkraken', 'length', 'info')
kraken=read_tsv(krakenfile, col_names=kraken_cols)[,2:3]

fullinfo.inclusive.kraken=full_join(fullinfo.inclusive, kraken, by='tig') %>%
    na.omit
speciesfile.inclusive.kraken=file.path(datadir, 'tigbin_species_inclusive_kraken.tsv')
write_tsv(fullinfo.inclusive.kraken, speciesfile.inclusive.kraken)


fullinfo.kraken=full_join(fullinfo, kraken, by='tig') %>%
    na.omit
speciesfile.kraken=file.path(datadir, 'tigbin_species_kraken.tsv')
write_tsv(fullinfo.kraken, speciesfile.kraken)
