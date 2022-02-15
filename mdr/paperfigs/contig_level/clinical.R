library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

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


cluster_copy(cluster, 'findMethFreq')

methgrouped=meth %>%
    filter(sum(methnum+umethnum)>5) %>%
    ##filter(sum(methnum+umethnum)>1) %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methfreq=methgrouped %>%
    do(findMethFreq(.))  %>%
    collect() %>%
    summarise(freq=mean(methfrac))

##keep only contigs that have every motif represented
nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])
methchroms=methfreq %>%
    rowwise() %>%
    filter(chrom %in% keepchroms)

##make matrix of meth info by contig
chrominfo=methchroms %>%
    spread(key=motif, value=freq)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom

##make dstance matrix
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

binnames=c()
for (i in rownames(matchrominfo)) {
    bin=chrombins$bin[chrombins$rname==i]
    binnames=c(binnames, bin)
}
numcolors=length(names(table(binnames)))-1
mycolors=c(colorRampPalette(brewer.pal(8, 'Set2'))(numcolors), '#000000')
colorkey=tibble(contig=rownames(table(binnames)), color=mycolors)
bincolors=c()
for (i in binnames) {
    color=colorkey$color[colorkey$contig==i]
    bincolors=c(bincolors, color)
}
bincolorinfo=data.frame(color=bincolors)
rownames(bincolorinfo)=rownames(matchrominfo)

##try dendrogram starting at 5.4 from tutorial below
##http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#ggplot2-integration
dend=matchrominfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
label_order=labels(dend)
label_colors=bincolorinfo[label_order,]
dend=dend %>%
    set('labels_col', label_colors)

chrombinspdf=file.path(dbxdir, 'clinical_contig_clusters_bin_colored.pdf')
pdf(chrombinspdf, h=8, w=30)
plot(dend)
dev.off()

binnedinfo=matchrominfo[binnames!='unknown',]
binneddend=binnedinfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
label_order=labels(binneddend)
label_colors=bincolorinfo[label_order,]
binneddend=binneddend %>%
    set('labels_col', label_colors)

chrombinspdf=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned.pdf')
pdf(chrombinspdf, h=8, w=30)
plot(binneddend)
dev.off()





####species info
tiginfocsv=file.path(dbxdir, 'tigbins_species.csv')
tiginfo=read_csv(tiginfocsv)


####attach classification info to the dendrogram
labelfull=NULL
for (i in label_order) {
    info=tiginfo[tiginfo$tig==i,]
    if (dim(info)[1]>0) {
        labelfull=c(labelfull, paste0(info$tig, ', ', info$tigleaf, ': ', info$bin, ', ',  info$binleaf, ' (', info$lcalevel, ')'))
    } else {
        ##if the contig wasn't assigned a taxononomy, but was assigned a bin
        bin=chrombins$bin[chrombins$rname==i]
        infobin=BAT[BAT$tig==paste0(bin, '.fasta'),][-(1:5)] %>% slice(1) %>% unlist(., use.names=FALSE)
        binindex=sum(!infobin=='no support')
        ##seems that one of the bins (bin_31) isn't even classified at the superkingdom level
        if (binindex>0) {
            binleaf=strsplit(infobin[binindex], ':', fixed=TRUE)[[1]][1]
        }else{
            binleaf='Organism'
        }
        labelfull=c(labelfull, paste0(i, ', not assigned: ', bin, ', ', binleaf, ' (NA)'))
        
    }
}

labels(binneddend)=labelfull

chrombinspdf=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull.pdf')
pdf(chrombinspdf, h=11, w=30)
par(mar = c(26, 4, 4, 2) + 0.1)
plot(binneddend)
dev.off()

##chrombinsvertpdf=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull_vert_test.pdf')
chrombinsvertpdf=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull_vert.pdf')
##pdf(chrombinsvertpdf, h=40, w=11)
pdf(chrombinsvertpdf, h=30, w=11)
par(mar = c(5, 4, 2, 40) + 0.1)
plot(binneddend, horiz=TRUE)
dev.off()

##large version
chrombinsvertlargepdf=file.path(dbxdir, 'clinical_contig_clusters_bin_colored_binned_labelfull_vert_large.pdf')
##pdf(chrombinsvertpdf, h=40, w=11)
pdf(chrombinsvertlargepdf, h=30, w=21)
par(mar = c(5, 4, 2, 40) + 0.1)
plot(binneddend, horiz=TRUE)
dev.off()



##check out weird contigs
treechrombins=chrombins %>%
    filter(rname %in% label_order)







####coverage hist analysis
covfile=file.path(projdir, 'mdr', 'align', paste0(prefix, 'polished.sorted.cov'))
cov_cols=c('tig', 'pos', 'cov')
cov=read_tsv(covfile, col_names=cov_cols) %>%
    group_by(tig) %>%
    mutate(covfrac=cov/max(cov))
cov_labeled=cov %>%
    filter(tig %in% label_order[1:20])

##get coverage spectrum
covfreq=cov %>%
    group_by(tig, cov) %>%
    summarise(freq=n()) %>%
    mutate(normfreq=freq/max(freq)) %>%
    mutate(normcov=cov/max(cov))
covfreq_labeled=covfreq %>%
    filter(tig %in% label_order)


covpeaks=covfreq %>%
    do(findpeaks(.)) %>%
    filter(tig %in% label_order) %>%
    mutate(composite=sum(h*t*s)) %>%
    arrange(-composite)


##obvs just do q-q plot??
covqq=covfreq %>%
    do(qqmse(.)) %>%
    filter(tig %in% label_order) %>%
    arrange(-mse)
##doesn't look as compeling as my weird thing looking for bimodality 
##tails are messing with me?

##plot hists
plotcov_head=covfreq %>%
    filter(tig %in% covpeaks_labeled$tig[1:20])
plotcov_tail=covfreq %>%
    filter(tig %in% tail(covpeaks_labeled$tig, 20))

histfile=file.path(dbxdir, 'clinical_contig_coverage_hists.pdf')
pdf(histfile, h=13, w=20)
headplot=ggplot(plotcov_head, aes(x=normcov, y=normfreq)) +
    geom_bar(stat='identity') +
    facet_wrap(~tig) +
    theme_bw()
print(headplot)
tailplot=ggplot(plotcov_tail, aes(x=normcov, y=normfreq)) +
    geom_bar(stat='identity') +
    facet_wrap(~tig) +
    theme_bw()
print(tailplot)
dev.off()





####plot clusters
clusters=cutree(binneddend, h=4)
clustinfo=tibble(tig=sapply(strsplit(names(clusters), ',', fixed=TRUE), '[[', 1),
                 cluster=paste0('cluster_',as.character(clusters))) %>%
    rowwise() %>%
    mutate(bin=chrombins$bin[chrombins$rname==tig]) %>%
    mutate(tiglen=chrombins$rlen[chrombins$rname==tig]) %>%
    rowwise() %>%
    mutate(binnum=as.numeric(strsplit(bin, '_', fixed=TRUE)[[1]][2])) %>%
    mutate(tigname=paste0(chartr('0123456789', 'abcdefghij', binnum), tig))

clustercolors=mycolors[1:15]
names(clustercolors)=names(table(clustinfo$bin))

plots=NULL
num=1
for (i in names(table(clustinfo$cluster))) {
    plotclustinfo=clustinfo %>%
        filter(cluster==i)
    plot=ggplot(plotclustinfo, aes(y=tig, x=tiglen, colour=bin, fill=bin)) +
        geom_bar(stat='identity', position='dodge', aes(alpha=.5)) +
        xlim(0, max(clustinfo$tiglen)) +
        ggtitle(i) +
        scale_fill_manual(values=clustercolors) +
        scale_color_manual(values=clustercolors) +
        scale_alpha(guide = 'none') +
        theme_bw()
    plots[[num]]=plot
    num=num+1
}

library(cowplot)    
clustplotfile=file.path(dbxdir, 'clinical_clusters.pdf')
pdf(clustplotfile, h=9, w=25)
##print(plot_grid(plotlist=plots, ncol=2))
plot=ggplot(clustinfo, aes(x=tigname, y=tiglen, colour=bin, fill=bin)) +
    geom_bar(stat='identity', position='dodge', alpha=.8, width=.85) +
    scale_fill_manual(values=clustercolors) +
    scale_color_manual(values=clustercolors) +
    facet_grid(~cluster, scales="free", space='free', switch='x') +
    scale_alpha(guide = 'none') +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.spacing=unit(1.5, 'lines'),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background=element_blank(),
          strip.placement='inside')
print(plot)
dev.off()




