library(tidyverse)
library(multidplyr)

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_asm'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'

methfile=file.path(datadir, 'clin_barocdes_methcalls.csv')
methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
meth=read_csv(methfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))


findMethFreq <- function(x) {
    ##collapses called meth into methfreq.
    ##x is meth df above, grouped by chrom and motif

    motiflen=nchar(x$motif[1])
    
    motifpos=NULL
    for (i in x$pos) {
        diffs=abs(x$pos-i)
        motifgroup=x[diffs<=motiflen,] %>%
            arrange(-methfrac) %>%
            mutate(totcalls=methnum+umethnum) %>%
            arrange(-totcalls)
        motifpos=bind_rows(motifpos, motifgroup[1,])
    }

    motifpos=unique(motifpos) %>%
        select(-totcalls)
    return(motifpos)
}
cluster_copy(cluster, 'findMethFreq')

methgrouped=meth %>%
    filter(sum(methnum+umethnum)>5) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methfreq=methgrouped %>%
    do(findMethFreq(.))  %>%
    collect() %>%
    summarise(freq=mean(methfrac))

nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])

methchroms=methfreq %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)

chrominfo=methchroms %>%
    spread(key=motif, value=freq)
