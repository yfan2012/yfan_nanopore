library(tidyverse)

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


##get idea of which tigs have which mtoifs represented
motifcounts=filter %>%
    group_by(chrom, motif) %>%
    summarise(num=n())


####isolate bases that are supposed to be methylated
##don't bother with the small staph plasmids since they only ahave a couple of motif calls total
basemeth=NULL
chrs=names(table(filter$chrom))
usechrs=chrs[1:(length(chrs)-2)]
for (i in usechrs) {
    label=chrlabels$label[chrlabels$chr==i]

    speciesmeth=filter %>%
        filter(chrom==i)

    for (j in 1:dim(methinfo)[1]) {
        qmotif=methinfo$motif[j]
        motifpos=methinfo$pos[j]
        motiflen=nchar(qmotif)

        
        speciesmotif=speciesmeth %>%
            filter(motif==qmotif)

        if (dim(speciesmotif)[1]>0) {
            basepos=seq(motifpos, dim(speciesmotif)[1], motiflen)
            methmotif=speciesmotif[basepos,] %>%
                mutate(label=paste0(motif, as.character(motifpos)))
            basemeth=bind_rows(basemeth, methmotif)
        }
    }
}

##plot
labels=names(table(basemeth$label))
basemethpdf=file.path(dbxdir, 'filter_motif_calls.pdf')
pdf(basemethpdf, h=9, w=15)
for (i in labels) {
    plotmeth=basemeth %>%
        filter(label==i)
    plot=ggplot(plotmeth, aes(x=chrom, y=methfrac, colour=chrom, fill=chrom, alpha=.5)) +
        geom_boxplot() +
        ggtitle(i) +
        theme_bw()
    print(plot)
}
dev.off()


####take top methfrac for each motif - does not assume any particular base is methylated
##very slow. parallel if running again
maxmeth=NULL
for (i in usechrs) {
    print(i)
    label=chrlabels$label[chrlabels$chr==i]

    speciesmeth=filter %>%
        filter(chrom==i)

    ##since just doing max, get rid of the redundant GATC motif
    methinfo2=methinfo[-9,]
    for (j in 1:dim(methinfo2)[1]) {
        qmotif=methinfo2$motif[j]
        print(qmotif)
        motifpos=methinfo2$pos[j]
        motiflen=nchar(qmotif)

        speciesmotif=speciesmeth %>%
            filter(motif==qmotif)
        
        if (dim(speciesmotif)[1]>0) {
            speciesmax=NULL
            for (k in 1:(dim(speciesmotif)[1]/motiflen)) {
                start=(k-1)*motiflen
                end=k*motiflen
                occurinfo=speciesmotif[start:end,]
                ##filter out lower confidence calls
                meancalls=(sum(occurinfo$methnum)+sum(occurinfo$umethnum))/2/motiflen
                validinfo=occurinfo %>%
                    filter((methnum+umethnum)>meancalls*.9)
                speciesmax=bind_rows(speciesmax, validinfo[validinfo$methfrac==max(validinfo$methfrac),])
            }
            maxmeth=bind_rows(maxmeth, speciesmax)
            
        }
    }
}

maxmethcsv=file.path(datadir, paste0(prefix, '.curate_maxcalls.csv'))
write_csv(maxmeth, maxmethcsv)

##plot maxmeth
motifs=names(table(maxmeth$motif))
maxmethpdf=file.path(dbxdir, 'filter_motif_maxcalls.pdf')
pdf(maxmethpdf, h=9, w=15)
for (i in motifs) {
    plotmeth=maxmeth %>%
        filter(motif==i)
    plot=ggplot(plotmeth, aes(x=chrom, y=methfrac, colour=chrom, fill=chrom, alpha=.5)) +
        geom_boxplot() +
        ggtitle(i) +
        theme_bw()
    print(plot)
}
dev.off()



####distances
maxmeth=read_csv(maxmethcsv)

##staph plasmid doesn't have full complement of motifs
maxmethagg=maxmeth %>%
    group_by(chrom, motif) %>%
    summarise(meanmeth=mean(methfrac))

nostaph=maxmethagg %>%
    filter(chrom!='Staphylococcus_aureus_plasmid1') %>%
    spread(key=motif, value=meanmeth)
matnostaph=as.matrix(nostaph %>% select(-chrom))
rownames(matnostaph)=nostaph$chrom
dists=dist(matnostaph)

abbmotif=maxmethagg %>%
    rowwise() %>%
    filter(motif %in% c('CAGAG', 'CCWGG', 'CTKVAG', 'GATC')) %>%
    spread(key=motif, value=meanmeth)
matabbmotif=as.matrix(abbmotif %>% select(-chrom))
rownames(matabbmotif)=abbmotif$chrom
abbdists=as.matrix(dist(matabbmotif))

abbdistsdf=as_tibble(abbdists) %>%
    mutate(chroms=rownames(abbdists)) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist,2))

##plot distances
distplotspdf=file.path(dbxdir, paste0('filter_motif_heatmaps_maxcalls.pdf'))
pdf(distplotspdf, h=9, w=9)
plot=ggplot(abbdistsdf, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill = dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
dev.off()


##same for basemeth
basemethagg=basemeth %>%
    group_by(chrom, motif) %>%
    summarise(meanmeth=mean(methfrac))
baseabbmotif=basemethagg %>%
    rowwise() %>%
    filter(motif %in% c('CAGAG', 'CCWGG', 'CTKVAG', 'GATC')) %>%
    spread(key=motif, value=meanmeth)
matbaseabbmotif=as.matrix(baseabbmotif %>% select(-chrom))
rownames(matbaseabbmotif)=baseabbmotif$chrom
baseabbdists=as.matrix(dist(matbaseabbmotif))

baseabbdistsdf=as_tibble(baseabbdists) %>%
    mutate(chroms=rownames(baseabbdists)) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist,2))

basedistplotspdf=file.path(dbxdir, paste0('filter_motif_heatmaps.pdf'))
pdf(basedistplotspdf, h=9, w=9)
plot=ggplot(baseabbdistsdf, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill = dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
dev.off()



##try for cmeth
cmethlist=c('GATC', 'CMTCGAKG', 'CCWGG', 'GCCGGC')
cagg=basemeth %>%
    group_by(chrom, motif) %>%
    summarise(meanmeth=mean(methfrac))
cabb=cagg %>%
    rowwise() %>%
    filter(motif %in% cmethlist) %>%
    spread(key=motif, value=meanmeth) %>%
    filter(chrom != 'Staphylococcus_aureus_plasmid1')
matcabb=as.matrix(cabb %>% select(-chrom))
rownames(matcabb)=cabb$chrom
cabbdists=as.matrix(dist(matcabb))

cabbdistsdf=as_tibble(cabbdists) %>%
    mutate(chroms=rownames(cabbdists)) %>%
    gather(key=chroms2, value=dist, -chroms) %>%
    mutate(rounded=round(dist,2))

cdistplotspdf=file.path(dbxdir, paste0('filter_motif_heatmaps_cmeth.pdf'))
pdf(cdistplotspdf, h=9, w=9)
plot=ggplot(cabbdistsdf, aes(x=chroms, y=chroms2)) +
    geom_tile(aes(fill = dist)) +
    geom_text(aes(label = rounded)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(plot)
dev.off()
