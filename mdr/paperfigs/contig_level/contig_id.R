library(tidyverse)

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'


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
    group_by(rname) %>%
    filter(n()>1) %>%
    summarise(all_bins=paste0(bin, collapse=','), all_rcov=paste0(total_rcov, collapse=','))

mumkeyfile=file.path(datadir, 'tigs2bins.tsv')
write_tsv(mumkey, mumkeyfile)
multikeyfile=file.path(datadir, 'tigs2bins_multi.tsv')
write_tsv(mummultiple, multikeyfile)



####look at organism classifications
CATfile=file.path(projdir, 'mdr/contig_id/CAT/200708_mdr_stool16native.CAT.names_official.txt')
phyloranks=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
CATcols=c('tig', 'classification', 'reason', 'lineage', 'lineage_scores', phyloranks)
CAT=read_tsv(CATfile)
names(CAT)=CATcols

BATfile=file.path(projdir, 'mdr/hiC/bin_id/BAT_single/200708_mdr_stool16native.BAT.names_official.txt')
BAT=read_tsv(BATfile)
names(BAT)=CATcols

tiginfo=NULL
for (i in CAT$tig[CAT$tig %in% mumkey$rname]) {
    tig=i
    bin=mumkey$bin[mumkey$rname==tig]

    
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

    supported=min(binindex, tigindex)
    lcamatches=which(allinfo$bin[1:supported]==allinfo$tig[1:supported])

    if (length(lcamatches)!=0) {
        lcaindex=max(lcamatches)
        lcalevel=phyloranks[lcaindex]
        lca=allinfo$tig[lcaindex]
    } else {
        lcalevel='none'
        lca='none'
    }
    
    tiginfo=bind_rows(tiginfo, tibble(tig=tig, bin=bin, tigleaf=tigleaf, binleaf=binleaf, lcalevel=lcalevel, lca=lca))
}

tiginfocsv=file.path(dbxdir, 'tigbins_species.csv')
write_csv(tiginfo, tiginfocsv)


binkey=unique(tiginfo %>% select(bin, binleaf))
