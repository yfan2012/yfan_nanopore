require(ggplot2)
require(Biostrings)
require(GenomicAlignments)
require(tidyverse)
require(gridExtra)
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

samps=c('samp15at24hrs','samp15at48hrs', 'samp20at24hrs', 'samp20at48hrs')
datadir='/kyber/Data/Nanopore/phage/'
##plot virus primary alignment lengths

alignall=foreach(i=samps, .combine=rbind) %dopar% {
    require(GenomicAlignments)
    bamfile=file.path(paste0(datadir,'align/',i, '.sorted.bam'))
    align=bam2df(bamfile)
    align$samp=i
    return(align)
}


    
pdf(file.path('~/Dropbox/Timplab_Data/phage/align/phage_align_lengths.pdf'), width=11, height=8.5)
print(ggplot(alignall, aes(x=samp, y=rlength, fill=samp)) +
      geom_violin(scale="width") +
      ggtitle('Phage Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
print(ggplot(alignall, aes(x=samp, y=rlength, fill=samp)) +
      geom_boxplot() +
      ggtitle('Phage Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
dev.off()
