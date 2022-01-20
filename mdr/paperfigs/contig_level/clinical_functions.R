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



findpeaks <- function(test) {
    ##normalized coverage frequency tibble
    ##return height and time (highest point of bimodality, and time spent in bimodality)
    
    steps=rev(seq(.01,1,.01))
    covrange=diff(range(test$cov))
    
    height=0
    time=0
    spread=0
    for (i in steps) {
        over=which(test$normfreq>=i)
        if (length(over)>1) {
            peaks=1+sum(diff(over)>1)
            if (peaks>1) {
                time=time+.01
                starts=c(over[1], over[which(diff(over)>1)]+1)
                ends=c(over[which(diff(over)>1)+1], over[length(over)])
                mids=starts+(ends-starts)/2
                ispread=max(diff(mids))/covrange
                if (spread<ispread && ispread>.05) {
                    spread=ispread
                    if (height<i){
                        height=i
                    }
                }
                
            }
        }
    }
    return(tibble(tig=test$tig[1], h=height, t=time, s=spread))
}


qqmse <- function(test){
    ##poisson coverage 
    mean=sum(test$cov*test$freq)/sum(test$freq)
    std=sqrt(mean)
    
    ordered_freq=sort(test$freq)
    ptile=test %>%
        rowwise() %>%
        mutate(percentile=sum(test$freq[test$cov<=cov])/sum(test$freq)) %>%
        mutate(theoretical=pnorm(cov, mean, std)) %>%
        mutate(res=(percentile=theoretical)^2)

    mse=sqrt(sum(ptile$res))
    return(tibble(tig=test$tig[1], mse=mse))
}


clustertigs <- function(methfreq, treefile, clusterfile) {
    
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

    ##add in bin info from mummer
    chrombinsfile='/mithril/Data/Nanopore/projects/methbin/paperfigs/contig_level/tigs2bins.tsv'
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
    
    ##look organism classifications
    dbxdir="~/gdrive/mdr/paperfigs/contig_level"
    tiginfocsv=file.path(dbxdir, 'tigbins_species.csv')
    tiginfo=read_csv(tiginfocsv)
    alltiginfo=full_join(tiginfo, tibble(tig=names(table(meth$chrom))), by='tig')
    alltiginfo[is.na(alltiginfo)]='unknown'
    
    ##attach classification info to the dendrogram
    projdir="/mithril/Data/Nanopore/projects/methbin"
    phyloranks=c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    CATcols=c('tig', 'classification', 'reason', 'lineage', 'lineage_scores', phyloranks)
    BATfile=file.path(projdir, 'mdr/hiC/bin_id/BAT_single/200708_mdr_stool16native.BAT.names_official.txt')
    BAT=read_tsv(BATfile)
    names(BAT)=CATcols
    
    labelfull=NULL
    for (i in label_order) {
        info=tiginfo[tiginfo$tig==i,]
        if (dim(info)[1]>0) {
            labelfull=c(labelfull, paste0(info$tig, ', ', info$tigleaf, ': ', info$bin, ', ',  info$binleaf, ' (', info$lca, ')'))
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
    
    
    pdf(treefile, h=30, w=11)
    par(mar = c(5, 4, 2, 40) + 0.1)
    plot(binneddend, horiz=TRUE)
    dev.off()
    
    ##plot clusters
    clusters=cutree(binneddend, h=5.2)
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
    

    pdf(clusterfile, h=9, w=25)
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
}
