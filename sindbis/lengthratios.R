require(ggplot2)
require(Biostrings)
require(GenomicAlignments)
require(tidyverse)
require(reshape2)
require(foreach)
require(doParallel)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)

##devel from plot for timp ubiome
bam2df <- function(bam) {
    dat.raw=readGAlignments(bam, use.names=T)
    dat.gr=granges(dat.raw)
    chroms=seqlevels(dat.gr)
    perchrom.len=data.frame(chr=factor(as.character(seqnames(dat.gr)), levels=chroms), rlength=width(dat.gr))

    if (dim(perchrom.len)[1] < 1) {
        perchrom.len=rbind(perchrom.len, data.frame(chr='none', rlength=0))
    }
    return(perchrom.len)
}

outdir='~/Dropbox/timplab_data/sindbis/plots/'
##excluding mock becuase it shouldn't have anything aligned
samps1=c('antibody', 'infected')
samps2=c('Antibody_2dpi', 'Antibody_3dpi', 'Sindbis_2dpi', 'Sindbis_3dpi')
allsamps=c(samps1, samps2)

##lengths to get ratios for
lens=seq(0,12000,100)

##get the ratios matrix for all samples, write to csv and plot to heatmap
for (samp in allsamps) {
    ##get alignment lengths
    bamfile=file.path(paste0('/dilithium/Data/Nanopore/sindbis/', samp,'/align/',samp,'.primary.sorted.bam'))
    aln=bam2df(bamfile)
    aln$samp=samp

    ##find how many reads are of each length
    numlens=foreach(i=1:length(lens), .combine=c) %dopar% {
        return(sum((aln$rlength<lens[i]+200) & (aln$rlength>lens[i]-200)))
    }

    ##get number of reads of each length with no overlap
    exactlens=foreach(i=1:length(lens), .combine=c) %dopar% {
        return(sum((aln$rlength<lens[i]+50) & (aln$rlength>lens[i]-50)))
    }
    
    fracs=data.frame(lengths=lens, numreads=numlens, frac=numlens/sum(numlens))
    exactfac=data.frame(lengths=lens, numreads=exactlens, frac=exactlens/sum(exactlens)) 
    write.table(fracs, paste0(outdir, samp, '_total_ratios.csv'), sep=',')
    write.table(exactfac, paste0(outdir, samp, '_total_ratios_nooverlap.csv'), sep=',')
    
    ##divide each length by ever other length
    ratios=foreach(i=1:length(numlens), .combine=rbind) %dopar% {
        return(numlens/numlens[i])
    }
    rownames(ratios)=as.character(lens)
    colnames(ratios)=as.character(lens)
    
    write.table(ratios, paste0(outdir, samp, '_length_ratios.csv'), row.names=TRUE, col.names=NA, sep=',')
    
    ratios[which(!is.finite(ratios))]=0
    ratiosdf=as.data.frame(ratios)
    ratiosdf$ypos=rownames(ratios)
    plotdf=melt(ratiosdf, id='ypos', value.name='heat')
    plotdf$ypos=as.numeric(plotdf$ypos)
    plotdf$variable=as.numeric(plotdf$variable)*100
    plotdf$log=abs(log(plotdf$heat))
    pdf(paste0(outdir, samp, '_length_ratios.pdf'), height=8.5, width=11)
    print(ggplot(plotdf, aes(x=variable, y=ypos)) +
        geom_tile(aes(fill=log)) +
        scale_fill_gradient2(low='red', high='blue') + 
        theme_bw())
    dev.off()
}
