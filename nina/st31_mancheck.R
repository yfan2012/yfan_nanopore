library(tidyverse)
library(RColorBrewer)

datadir='/pym/Data/Nanopore/projects/prolificans/st31/mancheck'
dbxdir='~/Dropbox/timplab_data/prolificans/mancheck'

pafnames=c('qname', 'qlen', 'qstart', 'qend', 'strand', 'rname', 'rlen', 'rstart', 'rend', 'matches', 'alen', 'mapq')
sca19tsv=file.path(datadir, 'st31.ragtag_fc.scaffolds.scaffold_19_RagTag.paf')
sca19=read_tsv(sca19tsv, col_names=pafnames)
tig22tsv=file.path(datadir, 'st31.ragtag_fc.scaffolds.contig_22_RagTag.paf')
tig22=read_tsv(tig22tsv, col_names=pafnames)

sca19common=sca19 %>%
    rowwise() %>%
    filter(qname %in% tig22$qname)
tig22common=tig22 %>%
    rowwise() %>%
    filter(qname %in% sca19$qname)
common=bind_rows(sca19common, tig22common)

##weight reads so that multiple alignments don't count as much
##this really just tells me that the ends in question are similar
weighted=common %>%
    group_by(qname, rname) %>%
    mutate(weight=1/n()) %>%
    rowwise() %>%
    mutate(pos=mean(rstart, rend)) %>%
    mutate(fracpos=pos/rlen)

st31_mergecheckfile=file.path(dbxdir, 'st31_merge.pdf')
pdf(st31_mergecheckfile, h=9, w=16)
ggplot(weighted, aes(x=fracpos, y=..density.., weight=weight, colour=rname, fill=rname, alpha=.3)) +
    geom_histogram(data=subset(weighted, rname=='scaffold_19_RagTag'), fill=brewer.pal(3, 'Set2')[1], colour=brewer.pal(3, 'Set2')[1], alpha=.3) +
    geom_histogram(data=subset(weighted, rname=='contig_22_RagTag'), fill=brewer.pal(3, 'Set2')[2], colour=brewer.pal(3, 'Set2')[2], alpha=.3) +
    ggtitle('Common reads') +
    theme_bw()
dev.off()


####check if any reads span
##get midpoints of alignments
span=common %>%
    filter(mapq!=0) %>% 
    group_by(qname, rname) %>%
    summarise(qpos=(qstart+qend)/2, len=mean(qlen)) %>%
    mutate(qfrac=qpos/len) ##%>%
    ##filter(qfrac>.55 || qfrac<.45)
uniquespan=unique(span) %>%
    group_by(qname) %>%
    filter(length(unique(rname))>1) %>% ##reads must align to both tigs
    filter(n()==2) %>% ##reads must have only two alignments
    filter(qfrac[1]-qfrac[2]>.4) ##midpoints of alignments must be 40% read length away
spanreads=unique(uniquespan$qname)

spaninfo=common %>%
    mutate(rpos=(rstart+rend)/2, rfrac=rpos/rlen) %>%
    rowwise() %>%
    filter(qname %in% spanreads) %>%
    filter(mapq==60) %>%
    group_by(qname) %>%
    filter(n()>1) %>%
    filter(abs(rfrac[1]-rfrac[2])>.8) %>%
    arrange(qname)




##look at ends
library(Biostrings)
library(GenomicRanges)
library(BSgenome)

fafile=file.path('/pym/Data/Nanopore/projects/prolificans/st31/genomes_final/st31.ragtag_fc.final.fasta')
asm=readDNAStringSet(fafile)

startranges=GRanges(seqnames=c('scaffold_19_RagTag', 'contig_22_RagTag'), ranges=IRanges(start=1, end=100))
endrange1=GRanges(seqnames=c('scaffold_19_RagTag'), ranges=IRanges(start=width(asm['scaffold_19_RagTag'])-100, end=width(asm['scaffold_19_RagTag'])))
endrange2=GRanges(seqnames=c('contig_22_RagTag'), ranges=IRanges(start=width(asm['contig_22_RagTag'])-100, end=width(asm['contig_22_RagTag'])))
starts=as.character(getSeq(asm, startranges))
end1=as.character(getSeq(asm, endrange1))
end2=as.character(getSeq(asm, endrange2))
