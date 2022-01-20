library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level'


####methylation distance
methfile=file.path(datadir, 'clin_barocdes_methcalls.perf.csv')
methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
meth=read_csv(methfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))

cluster_copy(cluster, 'findMethFreq')

methgrouped=meth %>%
    filter(sum(methnum+umethnum)>5) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methfreq=methgrouped %>%
    do(findMethFreq(.))  %>%
    collect() %>%
    summarise(freq=mean(methfrac))

treefile=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull_vert_perf.pdf')
clustplotfile=file.path(dbxdir, 'clinical_clusters_perf.pdf')
clustertigs(methfreq, treefile, clusterfile){




####coverage hist analysis
covfile=file.path(projdir, 'mdr', 'align', paste0(prefix, 'polished.perf.sorted.cov'))
cov_cols=c('tig', 'pos', 'cov')
cov=read_tsv(covfile, col_names=cov_cols) %>%
    group_by(tig) %>%
    mutate(covfrac=cov/max(cov))

tiglens=cov %>%
    group_by(tig) %>%
    summarise(len=max(pos))

##get coverage spectrum
covfreq=cov %>%
    group_by(tig, cov) %>%
    summarise(freq=n()) %>%
    mutate(normfreq=freq/max(freq)) %>%
    mutate(normcov=cov/max(cov))
    
covpeaks=covfreq %>%
    do(findpeaks(.)) %>%
    mutate(composite=sum(h*t*s)) %>%
    arrange(-composite) %>%
    mutate(len=tiglens$len[tiglens$tig==tig]) %>%
    ungroup() %>%
    mutate(rank=dense_rank(desc(composite)))



####motif frequency analysis
##try to rescue some contigs by eliminating non-useful motifs

##plot motif frequencies
motifdistpdf=file.path(dbxdir, 'meth_dists_per_motif.pdf')
pdf(motifdistpdf, h=10, w=17)
ggplot(methfreq, aes(x=freq, colour=motif, fill=motif, alpha=.5)) +
    geom_density() +
    facet_wrap(~motif, scales='free') +
    scale_color_manual(values=mycolors) +
    scale_fill_manual(values=mycolors) +
    theme_bw()
dev.off()

submethfreq=methfreq %>%
    filter(motif!='CTCGAG') %>%
    filter(motif!='GTCGAC')



####plasmid analysis
##make plasmid info
tigplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.plasmidfinder.tsv')
plas_cols=c('file', 'seq', 'start', 'end', 'strand', 'gene', 'coverage', 'covmap', 'gaps', 'covfrac', 'ident', 'db', 'acc', 'prod', 'res')
tigplas=read_tsv(tigplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc) %>%
    rowwise() %>%
    mutate(rank=covpeaks$rank[covpeaks$tig==seq]) %>%
    mutate(len=covpeaks$len[covpeaks$tig==seq]) %>%
    arrange(rank) %>%
    rename(tig=seq) %>%
    rename(tigident=ident) %>%
    mutate(mumbin=alltiginfo$bin[alltiginfo$tig==tig])


binplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.hiC.plasmidfinder.tsv')
binplas=read_tsv(binplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc) %>%
    rowwise() %>%
    mutate(seq=strsplit(seq, '.', fixed=TRUE)[[1]][1]) %>%
    rename(bin=seq) %>%
    rename(binident=ident)

plasinfo=full_join(tigplas, binplas, by='gene') %>%
    select(-tigident, -rank, -len, -binident) %>%
    arrange(gene)
##pause on this to include more plasmid tigs...



