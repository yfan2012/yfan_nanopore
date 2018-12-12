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
    require(GenomicAlignments)
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
datadir='/scratch/groups/mschatz1/cpowgs/phage/'
##plot sindbis primary alignment lengths

alignall=foreach(i=samps, .combine=rbind) %dopar% {
    require(GenomicAlignments)
    bamfile=file.path(paste0(datadir,'180809_phage/align/',i, '.sorted.bam'))
    align=bam2df(TEbam)
    return(align)
}

    
pdf(file.path('~/Dropbox/Timplab_Data/sindbis/sindbis_align_lengths.pdf'), width=11, height=8.5)
print(ggplot(allsamps, aes(x=samp, y=rlength, fill=samp)) +
      geom_violin(scale="width") +
      ggtitle('Sindbis Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
print(ggplot(allsamps, aes(x=samp, y=rlength, fill=samp)) +
      geom_boxplot() +
      ggtitle('Sindbis Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
dev.off()


##look at the spliced alignments to the genome to count the rat rna
ratab=file.path('/dilithium/Data/Nanopore/sindbis/antibody/align/antibody.rat.splicealn.primary.sorted.bam')
ratmock=file.path('/dilithium/Data/Nanopore/sindbis/mock/align/mock.rat.splicealn.primary.sorted.bam')
ratTE=file.path('/dilithium/Data/Nanopore/sindbis/infected/align/infected.rat.splicealn.primary.sorted.bam')

ratabcount=bam2df(ratab)
ratmockcount=bam2df(ratmock)
ratTEcount=bam2df(ratTE)

