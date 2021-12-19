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

##keep only chroms that have every motif represented
nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])
methchroms=methfreq %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)

chrominfo=methchroms %>%
    spread(key=motif, value=freq)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom

chromdists=as_tibble(as.matrix(dist(matchrominfo))) %>%
    mutate(chroms=chrominfo$chrom) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist, 2))

chromdistspdf=file.path(dbxdir, 'clinical_contig_distances.pdf')
pdf(chromdistspdf, h=55, w=55)
plot=ggplot(chromdists, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill=dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
dev.off()

chromclustspdf=file.path(dbxdir, 'clinical_contig_clusters.pdf')
pdf(chromclustspdf, h=8, w=30)
plot(hclust(dist(matchrominfo)))
print(plot)
dev.off()



####add in bin info from mummer
chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)
for (i in keepchroms) {
    if (!(i %in% chrombins$rname)) {
        info=tibble(rname=i, bin='unknown', total_rcov=NA)
        chrombins=bind_rows(chrombins, info)
    }
}
methbins=methchroms %>%
    rowwise() %>%
    mutate(bin=chrombins$bin[chrombins$rname==chrom])
