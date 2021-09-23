library(tidyverse)
library(Biostrings)
source('~/Code/yfan_nanopore/mdr/disco/contig_id_functions.R')

prefix='200708_mdr_stool16native_polished'
projdir='/mithril/Data/Nanopore/projects/methbin/mdr'


####plasmidfinder
plasmidtsv=file.path(projdir, 'amr', paste0('200708_mdr_stool16native.plasmidfinder.tsv'))
pf_cols=c('file', 'tig', 'start', 'end', 'strand', 'gene', 'covinfo', 'covmap', 'gaps', 'cov', 'ident', 'db', 'acc', 'resistance')
plasmid=read_tsv(plasmidtsv, col_names=pf_cols, comment='#') %>%
    select(c(-file, -strand, -covinfo, -covmap, -gaps, -db, -resistance)) %>%
    filter(ident>99) %>%
    mutate(overall=cov*ident) %>%
    group_by(tig) %>%
    mutate(max=max(overall)) %>%
    ungroup() %>%
    filter(overall==max) %>%
    select(c(-max, -overall))


####blast db info
faheaderfile=gzfile('/atium/Data/ref/bacteria/all_bacteria_refs_faheaders.txt')
accinfo=read_table2(faheaderfile, col_names=c('acc', 'genus', 'species'))



####blast against bacteria ref
blasttsv=file.path(projdir, 'blast_contigs', paste0(prefix, '.assembly.tsv'))
blastcols=c('tig', 'ref', 'ident', 'alen', 'mismatch', 'gaps', 'tigstart', 'tigend', 'refstart', 'refend', 'eval', 'bit')
blast=read_tsv(blasttsv, col_names=blastcols, comment='#') %>%
    select(-c(mismatch, gaps, refstart, refend, eval, bit)) %>%
    filter(alen>10000 & ident>98) %>%
    filter(!tig %in% plasmid$tig)
merged=blast %>%
    group_by(tig, ref) %>%
    do(merge_overlaps(.))

mergecollapse=merged %>%
    group_by(tig, ref) %>%
    summarise(totalbases=sum(end-start))

tigidinfo=mergecollapse %>%
    rowwise() %>%
    mutate(place=which(ref==accinfo$acc)) %>%
    mutate(taxid=accinfo$taxid[place]) %>%
    rowwise() %>%
    mutate(place=which(accinfo$acc==ref)) %>%
    mutate(g=accinfo$genus[place], s=accinfo$species[place]) %>%
    select(-place)
tigidinfofile=file.path(projdir, 'contig_id', 'blast_bacter_ref.csv')
write_csv(tigidinfo, tigidinfofile)


tigid=tigidinfo %>%
    group_by(tig) %>%
    summarise(refmax=ref[which(totalbases==max(totalbases))],
              basemax=totalbases[which(totalbases==max(totalbases))],
              gmax=g[which(totalbases==max(totalbases))],
              smax=s[which(totalbases==max(totalbases))]) %>%
    arrange(-basemax)
tigidfile=file.path(projdir, 'contig_id', 'blast_tigid.csv')
write_csv(tigid, tigidfile)
