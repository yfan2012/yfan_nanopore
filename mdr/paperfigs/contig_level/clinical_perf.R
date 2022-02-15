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
##methfile=file.path(datadir, 'clin_barocdes_methcalls.perf.csv')
methfile=file.path(datadir, 'clin_barocdes_methcalls.perf2.csv')
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


##try to rescue some contigs by eliminating non-useful motifs
##plot motif frequencies
motifdistpdf=file.path(dbxdir, 'meth_dists_per_motif_bc2.pdf')
##motifdistpdf=file.path(dbxdir, 'meth_dists_per_motif.pdf')
pdf(motifdistpdf, h=10, w=17)
mycolors=colorRampPalette(brewer.pal(8, 'Set2'))(length(table(methfreq$motif)))
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

treefile=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull_vert_perf_13motifs.pdf')
clusterfile=file.path(dbxdir, 'clinical_clusters_perf_13motifs.pdf')
clustertigs(submethfreq, treefile, clusterfile, 4)


exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
testmethfreq=methfreq %>%
    ##exlude based on discrim power from distributions above
    rowwise() %>%
    filter(!motif %in% exvec)

freqs=testmethfreq %>%
    spread(key=motif, value=freq)
nacount=colSums(is.na(freqs))/dim(freqs)[1]

treefile=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull_vert_perf_test.pdf')
clusterfile=file.path(dbxdir, 'clinical_clusters_perf_test.pdf')
heatfile=file.path(dbxdir, 'clinical_clusters_heatmap_perf_test.pdf')
clustertigs(testmethfreq, treefile, clusterfile, 4)
heatcheck(freqs, heatfile)




####coverage hist analysis
covfile=file.path(projdir, 'mdr', 'align', paste0('200708_mdr_stool16native_asmpolished.perf.sorted.cov'))
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







####plasmid analysis
##make plasmid info
chrombinsfile=file.path(datadir, 'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

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
    mutate(mumbin='unknown') %>%
    mutate(genebin='unknown')
for (i in tigplas$tig) {
    if (i %in% chrombins$rname) {
        tigplas$mumbin[tigplas$tig==i]=chrombins$bin[chrombins$rname==i]
    }
}

binplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.hiC.plasmidfinder.tsv')
binplas=read_tsv(binplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc) %>%
    rowwise() %>%
    mutate(seq=strsplit(seq, '.', fixed=TRUE)[[1]][1]) %>%
    rename(bin=seq) %>%
    rename(binident=ident)

for (i in tigplas$gene) {
    if (i %in% binplas$gene) {
        bins=paste(unique(binplas$bin[binplas$gene==i]), collapse='|')
        tigplas$genebin[tigplas$gene==i]=bins
    }
}

##get contigs nearest to plasmids
plasnear=NULL
for (i in tigplas$tig) {
    plasnearest=plascheck(freqs, i)
    plasnear=bind_rows(plasnear, plasnearest[1:20,])
}
plasnear=unique(plasnear[complete.cases(plasnear),]) %>%
    select(-rounded)

##see if nearest contigs have a bin assignment
plasnearbins=NULL
for (i in 1:dim(plasnear)[1]) {
    info=plasnear[i,]

    if (info$chroms2 %in% chrombins$rname) {
        info$bin=chrombins$bin[chrombins$rname==info$chroms2]
    } else {
        info$bin='unknown'
    }
    plasnearbins=bind_rows(plasnearbins, info)
}



####classification methods
##nearest only
nearest=plasnearbins %>%
    filter(bin!='unknown') %>%
    group_by(chroms) %>%
    filter(row_number()==1) %>%
    mutate(mumbin=paste(unique(tigplas$mumbin[tigplas$tig==chroms], collapse='|')))

##nearest known
nearestknown=plasnearbins %>%
    group_by(chroms) %>%
    do(plas_nearest_known(.,3)) %>%
    mutate(mumbin=paste(unique(tigplas$mumbin[tigplas$tig==chroms]), collapse='|')) %>%
    mutate(genebin=paste(unique(tigplas$genebin[tigplas$tig==chroms]), collapse='|')) %>%
    filter(nummotifs>=10) %>%
    mutate_if(is.character, str_replace_all, pattern='nothing1', replacement='none') %>%
    mutate_if(is.character, str_replace_all, pattern='nothing2', replacement='none') %>%
    mutate_if(is.character, str_replace_all, pattern='nothing3', replacement='none')

##add in taxonomy info
tiginfocsv=file.path(dbxdir, 'tigbins_species.csv')
tiginfo=read_csv(tiginfocsv)

plasneartax=plasnearbins %>%
    group_by(chroms, chroms2) %>%
    do(taxonomy_pairs(.))







####rand distance
##plain dendrogram labels to make this easier
plaindend=binnedinfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
##truth
truthbins=tibble(tig=labels(plaindend)) %>%
    rowwise() %>%
    filter(tig %in% tiginfo$tig) %>%
    mutate(bin=tiginfo$bin[which(tiginfo$tig==tig)])


clusts=cutree(plaindend, h=5)
