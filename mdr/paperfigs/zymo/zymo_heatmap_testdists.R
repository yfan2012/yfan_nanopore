library(tidyverse)
library(multidplyr)
library(pheatmap)
library(GenomicRanges)
library(Biostrings)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')
cluster_copy(cluster, 'findMethFreq')

projdir='/mithril/Data/Nanopore/projects/methbin/zymo'
prefix='20190809_zymo_control'
datadir=file.path(projdir, 'contig_agg', prefix)
filterfile=file.path(datadir, paste0(prefix,'.curate_filter.csv'))

dbxdir='~/gdrive/mdr/zymo'


##read motif info
cmethfile=file.path(projdir, 'truth/bisulfite/zymo_cmeth.csv')
cmeth=read_csv(cmethfile)
amethfile=file.path(projdir, 'truth/pacbio/zymo_ameth.csv')
ameth=read_csv(amethfile)
methinfo=bind_rows(cmeth,ameth) %>%
    group_by(motif, pos) %>%
    summarise(num=n())
    
nametochrfile='~/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
chrlabels=read_table2(nametochrfile, col_names=c('chr', 'label'))

##read filtered megalodon output
filter_cols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
filter=read_csv(filterfile, col_names=filter_cols) %>%
    filter(!grepl('tig0', chrom, fixed=TRUE)) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))

##count motif occurences themselves
reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
ref=readDNAStringSet(reffile, format='fasta')
refnames=sapply(strsplit(names(ref), ' '), function(x) x[1])
names(ref)=refnames

##get idea of which tigs have which mtoifs represented
motifcounts=filter %>%
    group_by(chrom, motif) %>%
    summarise(num=n())




####isolate bases that are supposed to be methylated
##don't bother with the small staph plasmids since they only ahave a couple of motif calls total
aggregate_meth  <- function(speciesmeth, methinfo){
    basemeth=NULL
    for (j in 1:dim(methinfo)[1]) {
        qmotif=methinfo$motif[j]
        print(qmotif)
        motifpos=methinfo$pos[j]
        motiflen=nchar(qmotif)
        
        speciesmotif=speciesmeth %>%
            filter(motif==qmotif)
        
        if (dim(speciesmotif)[1]>0) {
            speciesfreq=methFreqByPos(speciesmotif, motifpos) %>%
                mutate(label=paste0(motif, as.character(motifpos)))
            basemeth=bind_rows(basemeth, speciesfreq)
        }
    }
    return(basemeth)
}

cluster_copy(cluster, 'aggregate_meth')
cluster_copy(cluster, 'methinfo')
cluster_copy(cluster, 'methFreqByPos')

chrs=names(table(filter$chrom))
usechrs=chrs[1:(length(chrs)-2)]

speciesmeth=filter %>%
    filter(chrom %in% usechrs) %>%
    group_by(chrom) %>%
    partition(cluster)
basemeth=speciesmeth %>%
    do(aggregate_meth(., methinfo)) %>%
    collect() %>%
    ungroup()


##print out number of motif occurences and number of 
motiflist=names(table(basemeth$label))
chromslist=names(table(filter$chrom))
motifcounts=tibble(chrom=rep(chromslist, each=length(motiflist)),
                   motif=rep(motiflist, times=length(chromslist))) %>%
    mutate(motif=substr(motif, start=1, stop=str_length(motif)-1)) %>%
    rowwise() %>%
    mutate(total=countPattern(motif, ref[[chrom]], fixed=FALSE)) %>%
    distinct()
coveredcounts=basemeth %>%
    group_by(chrom, label) %>%
    summarise(num=n()) %>%
    mutate(motif=substr(label, start=1, stop=str_length(label)-1)) %>%
    select(-label) %>%
    distinct()
allcounts=full_join(motifcounts, coveredcounts, by=c('chrom', 'motif')) %>%
    mutate(frac=num/total)
allcounts[is.na(allcounts)]=0



##mess with dists
methagg=basemeth %>%
    group_by(chrom, label) %>%
    summarise(meanmeth=mean(methfrac))

abbmotif=methagg %>%
    rowwise() %>%
    filter(label %in% c('GATC2', 'CAGAG4', 'GATC4', 'CCWGG2', 'CTKVAG5')) %>%
    ##filter(label %in% c('GATC2', 'CAGAG4', 'GATC4', 'CCWGG2', 'CTKVAG5', 'CMTCGAKG4', 'GCCGGC3', 'TCCGGA3')) %>%
    spread(key=label, value=meanmeth)
matabbmotif=as.matrix(abbmotif %>% select(-chrom))
rownames(matabbmotif)=abbmotif$chrom
##abbdists=as.matrix(dist(matabbmotif)) ##all allowed motifs
##abbdists=as.matrix(dist(matabbmotif, method='manhattan'))
abbdists=as.matrix(dist(matabbmotif, method='canberra'))

abbdistsdf=as_tibble(abbdists) %>%
    mutate(chroms=rownames(abbdists)) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist,2))

clustfeats=as.data.frame(abbmotif[,2:6])
rownames(clustfeats)=abbmotif$chrom
chrorder=labels(as.dendrogram(hclust(dist(clustfeats), 'average')))
chrlabs=tibble(newlab=paste0(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'), chrorder),
               oldlab=chrorder)

abbdistsdf=abbdistsdf %>%
    rowwise() %>%
    mutate(chroms=chrlabs$newlab[chrlabs$oldlab==chroms]) %>%
    mutate(chroms2=chrlabs$newlab[chrlabs$oldlab==chroms2])


##merged motifs
motifmethagg=basemeth %>%
    group_by(chrom, motif) %>%
    summarise(meanmeth=mean(methfrac))
cleanmotif=motifmethagg %>%
    rowwise() %>%
    filter(motif %in% c('GATC', 'CAGAG', 'CCWGG')) %>%
    ##filter(motif %in% c('GATC', 'CAGAG', 'CCWGG', 'CMTCGAKG', 'GCCGGC', 'TCCGGA')) %>%
    spread(key=motif, value=meanmeth)
matcleanmotif=as.matrix(cleanmotif %>% select(-chrom))
rownames(matcleanmotif)=cleanmotif$chrom
##cleandists=as.matrix(dist(matcleanmotif))
##cleandists=as.matrix(dist(matcleanmotif, method='manhattan'))
cleandists=as.matrix(dist(matcleanmotif, method='canberra'))

cleandistsdf=as_tibble(cleandists) %>%
    mutate(chroms=rownames(cleandists)) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist,2))

chrlabs=tibble(newlab=paste0(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'), chrorder),
               oldlab=chrorder)

cleandistsdf=cleandistsdf %>%
    rowwise() %>%
    mutate(chroms=chrlabs$newlab[chrlabs$oldlab==chroms]) %>%
    mutate(chroms2=chrlabs$newlab[chrlabs$oldlab==chroms2])


##plot
##distplotspdf=file.path(dbxdir, paste0('heatmap_testdist.pdf'))
##distplotspdf=file.path(dbxdir, paste0('heatmap_testmanhat.pdf'))
distplotspdf=file.path(dbxdir, paste0('heatmap_testcanberra.pdf'))
pdf(distplotspdf, h=9, w=9)
plot=ggplot(abbdistsdf, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill = dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "#FC8D62", high = "white") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
plot=ggplot(cleandistsdf, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill = dist)) +
    geom_text(aes(label = rounded)) +
    ggtitle('Merged motifs') +
    scale_fill_gradient(low = "#FC8D62", high = "white") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
dev.off()



