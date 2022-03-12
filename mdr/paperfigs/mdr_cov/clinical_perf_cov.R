library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(16)
cluster_library(cluster, 'tidyverse')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/contig_level_cov'


####methylation distance
##methfile=file.path(datadir, 'clin_barocdes_methcalls.perf.csv')
methfile=file.path(datadir, 'clin_barocdes_methcalls.perf2.csv')
methcols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
meth=read_csv(methfile, col_names=methcols) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))

##coverage
covcols=c('chrom', 'pos', 'cov')
covfile=file.path(projdir, 'mdr/align/200708_mdr_stool16native_asmpolished.perf.sorted.cov')
cov=read_tsv(covfile, covcols)

methcov=inner_join(meth, cov, by=c('chrom', 'pos'))



cluster_copy(cluster, 'findMethFreq')

methgrouped=methcov %>%
    filter(sum(methnum+umethnum)>5) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methloci=methgrouped %>%
    do(findMethFreq(.))  %>%
    collect()
methfreq=methloci %>%
    ungroup() %>%
    mutate(methfrac=methnum/cov) %>%
    group_by(motif, chrom) %>%
    summarise(freq=mean(methfrac), numloci=n())



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
    select(-numloci) %>%
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
tiginfocsv=file.path('~/gdrive/mdr/paperfigs/contig_level/tigbins_species.csv')
tiginfo=read_csv(tiginfocsv)

plasneartax=plasnearbins %>%
    group_by(chroms, chroms2) %>%
    do(taxonomy_pairs(.))

nearestcsv=file.path(dbxdir, 'nearestplas_cov.csv')
write_csv(nearest, nearestcsv)
nearestknowncsv=file.path(dbxdir, 'nearestknown_cov.csv')
write_csv(nearestknown, nearestknowncsv)





####rand distance and contamination analysis
##from clustertigs function, but stops at the plain dend
methfreq=testmethfreq
nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])
methchroms=methfreq %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)
chrominfo=methchroms %>%
    spread(key=motif, value=freq)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom
plaindend=matchrominfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram

##truth
truthbins=tibble(tig=labels(plaindend)) %>%
    filter(tig %in% tiginfo$tig) %>%
    rowwise() %>%
    mutate(bin=chrombins$bin[which(chrombins$rname==tig)]) %>%
    mutate(tiglen=chrombins$rlen[chrombins$rname==tig]) %>%
    filter(bin!='unknown')

library(fossil)
library(mclust)
maxheight=attributes(plaindend)$height
rands=NULL
contams=NULL

for (height in seq(0, maxheight, .01)){
    clusts=cutree(plaindend, h=height)
    truthbins=truthbins %>%
        rowwise() %>%
        mutate(clustbin=clusts[tig])
    numclusts=length(table(truthbins$clustbin))
    
    ##rand index stuff
    ri=rand.index(as.numeric(as.factor(truthbins$bin)), truthbins$clustbin)
    ari=adjustedRandIndex(truthbins$bin, truthbins$clustbin)
    heightrand=tibble(height=height,
                      numclusts=numclusts,
                      ri=ri,
                      ari=ari)
    rands=bind_rows(rands, heightrand)

    ##togetherness - what percentage of contigs is not in the right cluster?
    ##assign bins to clusters by majority
    majclusts=truthbins %>%
        group_by(bin, clustbin) %>%
        summarise(seq=sum(tiglen)) %>%
        filter(seq==max(seq))
    ##for each contig, see if its part of the 'right' cluster
    contam=truthbins %>%
        rowwise() %>%
        mutate(inmajority=clustbin==majclusts$clustbin[majclusts$bin==bin])
    majority=sum(contam$tiglen[contam$inmajority])/sum(contam$tiglen)
    nummajority=sum(contam$inmajority)/dim(contam)[1]
    
    ##purity - what percentage of a meth cluster is not the majority bin of that cluster?
    ##assign clusters to bins according to majority composition
    pureclusts=truthbins %>%
        group_by(clustbin, bin) %>%
        summarise(seq=sum(tiglen)) %>%
        ungroup() %>%
        group_by(clustbin) %>%
        filter(seq==max(seq))
    ##for each contig, see if it's part of the 'right' cluster
    purity=truthbins %>%
        rowwise() %>%
        mutate(purity=bin==pureclusts$bin[pureclusts$clustbin==clustbin])
    pure=sum(purity$tiglen[purity$purity])/sum(purity$tiglen)
    numpure=sum(purity$purity)/dim(purity)[1]

    heightcontam=tibble(height=height,
                      numclusts=numclusts,
                      majority=majority,
                      revmajority=1-majority,
                      nummajority=nummajority,
                      revnummajority=1-nummajority,
                      pure=pure,
                      numpure=numpure)
    contams=bind_rows(contams, heightcontam)
}


rands=rands %>%
    select(-height) %>%
    unique() %>%
    arrange(numclusts)

ucontams=contams %>%
    select(-height) %>%
    unique() %>%
    arrange(numclusts)
percent=ucontams %>%
    select(-c(nummajority, numpure, revmajority, revnummajority)) %>%
    gather('key', 'value', -numclusts)
num=ucontams %>%
    select(-c(majority, pure, revmajority, revnummajority)) %>%
    gather('key', 'value', -numclusts)


metricpdf=file.path(dbxdir, 'clinical_contig_clusters_metrics_perf.pdf')
pdf(metricpdf, h=8, w=11)
percentplot=ggplot(percent, aes(x=numclusts, y=value, colour=key)) +
    geom_line() +
    ggtitle('Percent sequence') +
    theme_bw()
rocplot=ggplot(ucontams, aes(x=revmajority, y=pure)) +
    geom_step() +
    ggtitle('Percent seqeunce') +
    xlim(0,1) +
    ylim(0,1) +
    geom_abline(a=0, b=1) +
    theme_bw()
numplot=ggplot(num, aes(x=numclusts, y=value, colour=key)) +
    geom_line() +
    ggtitle('Percent contigs') +
    theme_bw()
rocnumplot=ggplot(ucontams, aes(x=revnummajority, y=numpure)) +
    geom_step() +
    ggtitle('Percent contigs') +
    xlim(0,1) +
    ylim(0,1) +
    geom_abline(a=0, b=1) +
    theme_bw()
plot(percentplot)
plot(rocplot)
plot(numplot)
plot(rocnumplot)
dev.off()

ratings=contams %>%
    mutate(seqrate=majority*pure) %>%
    mutate(tigrate=nummajority*numpure) %>%
    unique() %>%
    select(-c(height))

##check randomization under the same structure - randomize labels on the existing tree structure
truthbins=tibble(tig=labels(plaindend)) %>%
    rowwise() %>%
    filter(tig %in% tiginfo$tig) %>%
    mutate(bin=tiginfo$bin[which(tiginfo$tig==tig)]) %>%
    mutate(tiglen=chrombins$rlen[chrombins$rname==tig])

allroc=NULL
for (i in 1:50) {
    randchrominfo=as.matrix(chrominfo %>% select(-chrom))
    rownames(randchrominfo)=sample(chrominfo$chrom)
    randodend=randchrominfo %>%
        scale %>% 
        dist %>%
        hclust %>%
        as.dendrogram
    randroc=get_tree_roc(randodend, truthbins) %>%
        mutate(samp=i)
    allroc=bind_rows(allroc, randroc)
}

numtips=length(chrominfo$chrom)
allrands=NULL
for (i in 1:50) {
    randchrominfo=matrix(sample(matchrominfo), nrow=dim(matchrominfo)[1], ncol=dim(matchrominfo)[2])
    rownames(randchrominfo)=chrominfo$chrom
    colnames(randchrominfo)=colnames(matchrominfo)
    randodend=randchrominfo %>%
        scale %>% 
        dist %>%
        hclust %>%
        as.dendrogram
    randroc=get_tree_roc(randodend, truthbins) %>%
        mutate(samp=i)
    allrands=bind_rows(allrands, randroc)
}

realroc=get_tree_roc(plaindend, truthbins) %>%
    mutate(samp=0)

plotallroc=bind_rows(realroc, allroc) %>%
    mutate(label=case_when(samp==0 ~ 'real', TRUE ~ 'rando')) %>%
    mutate(seqtogether=1-seqtogether) %>%
    mutate(numtogether=1-numtogether)
plotallrands=bind_rows(realroc, allrands) %>%
    mutate(label=case_when(samp==0 ~ 'real', TRUE ~ 'rando')) %>%
    mutate(seqtogether=1-seqtogether) %>%
    mutate(numtogether=1-numtogether)

plotallroccsv=file.path(dbxdir, 'plotallroc_cov.csv')
write_csv(plotallroc, plotallroccsv)
plotallrandscsv=file.path(dbxdir, 'plotallrands_cov.csv')
write_csv(plotallrands, plotallrandscsv)


##randopdf=file.path(dbxdir, 'clinical_contig_clusters_metrics_perf_rando.pdf')
randopdf=file.path(dbxdir, 'clinical_contig_clusters_metrics_perf_rando_points.pdf')
pdf(randopdf, h=8, w=11)
plot1=ggplot(plotallroc %>% filter(label=='real'), aes(x=seqtogether, y=seqpure, colour=label, alpha=.02))+
##plot1=ggplot(plotallroc, aes(x=seqtogether, y=seqpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallroc %>% filter(label=='rando'), mapping=aes(x=seqtogether, y=seqpure, colour=label, alpha=.02)) +
    ggtitle('Percent Sequence: random labels, same tree structure') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot1)
plot2=ggplot(plotallrands %>% filter(label=='real'), aes(x=seqtogether, y=seqpure, colour=label, alpha=.02))+
##plot2=ggplot(plotallrands, aes(x=seqtogether, y=seqpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallrands %>% filter(label=='rando'), mapping=aes(x=seqtogether, y=seqpure, colour=label, alpha=.02)) +
    ggtitle('Percent Sequence: random structure') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot2)
plot3=ggplot(plotallroc %>% filter(label=='real'), aes(x=numtogether, y=numpure, colour=label, alpha=.02))+
##plot3=ggplot(plotallroc, aes(x=numtogether, y=numpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallroc %>% filter(label=='rando'), mapping=aes(x=numtogether, y=numpure, colour=label, alpha=.02)) +
    ggtitle('Percent Contigs: random labels, same tree structure') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot3)
plot4=ggplot(plotallrands %>% filter(label=='real'), aes(x=numtogether, y=numpure, colour=label, alpha=.02))+
##plot4=ggplot(plotallrands, aes(x=numtogether, y=numpure, colour=label, alpha=.02))+
    geom_step() +
    geom_point(plotallrands %>% filter(label=='rando'), mapping=aes(x=numtogether, y=numpure, colour=label, alpha=.02)) +
    ggtitle('Percent Contigs: random sturcture') +
    scale_colour_brewer(palette='Set2') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot4)
dev.off()




####investigate why cov gets more wrong
set1=c('contig_16', 'contig_1253', 'contig_386')
set2=c('contig_433', 'contig_16', 'contig_774')
set3=c('contig_503', 'contig_161', 'contig_137')

setloci=methloci %>%
    filter(chrom %in% set1)
covfreq=setloci %>%
    ungroup() %>%
    mutate(methfrac=methnum/cov) %>%
    group_by(motif, chrom) %>%
    summarise(freq=mean(methfrac), numloci=n()) %>%
    select(-numloci) %>%
    ##exlude based on discrim power from distributions above
    rowwise() %>%
    filter(!motif %in% exvec)
covmat=covfreq %>%
    spread(key=motif, value=freq)

maxfreq=setloci %>%
    ungroup() %>%
    group_by(motif, chrom) %>%
    summarise(freq=mean(methfrac), numloci=n()) %>%
    select(-numloci) %>%
    rowwise() %>%
    filter(!motif %in% exvec)
maxmat=maxfreq %>%
    spread(key=motif, value=freq)

squashmotifs=c('KCCGGM', 'GGCC', 'GCGC')
squashloci=setloci %>%
    filter(motif %in% squashmotifs) %>%
    mutate(methcov=methnum/cov) %>%
    mutate(umethcov=umethnum/cov)

squashpdf=file.path(dbxdir, 'clinical_contig_squash.pdf')
pdf(squashpdf, h=8, w=17)
for (i in squashmotifs) {
    squashmotif=squashloci %>%
        filter(motif==i)
    numpoints=min(table(squashmotif$chrom))
    plotsquash=squashmotif %>%
        group_by(chrom) %>%
        filter(row_number()<=numpoints)
    plot=ggplot(plotsquash, aes(x=umethcov, y=methcov, colour=chrom, alpha=.05)) +
        geom_point() +
        facet_wrap(~chrom) +
        ggtitle(i) +
        scale_colour_brewer(palette='Set2') +
        theme_bw()
    print(plot)
}
dev.off()
