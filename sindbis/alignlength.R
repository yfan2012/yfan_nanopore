require(ggplot2)
require(Biostrings)
require(GenomicAlignments)
require(tidyverse)
require(gridExtra)

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


abbam=file.path('/dilithium/Data/Nanopore/sindbis/antibody/align/antibody.sorted.bam')
mockbam=file.path('/dilithium/Data/Nanopore/sindbis/mock/align/mock.sorted.bam')
TEbam=file.path('/dilithium/Data/Nanopore/sindbis/infected/align/infected.sorted.bam')
ab=bam2df(abbam)
ab$samp='antibody'
mock=bam2df(mockbam)
mock$samp='mock'
TE=bam2df(TEbam)
TE$samp='infected'

allsamps=rbind(ab, mock, TE)


pdf(file.path('~/Dropbox/Timplab_Data/sindbis/align_lengths.pdf'), width=11, height=8.5)
print(ggplot(allsamps, aes(x=samp, y=rlength, fill=samp)) +
      geom_violin(scale="width") +
      ggtitle('Alignment lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
print(ggplot(allsamps, aes(x=samp, y=rlength, fill=samp)) +
      geom_boxplot() +
      ggtitle('Alignment lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
dev.off()


ratab=file.path('/dilithium/Data/Nanopore/sindbis/antibody/align/antibody.rat.splicealn.sorted.bam')
ratmock=file.path('/dilithium/Data/Nanopore/sindbis/mock/align/mock.rat.splicealn.sorted.bam')
ratTE=file.path('/dilithium/Data/Nanopore/sindbis/infected/align/infected.rat.splicealn.sorted.bam')

ratabcount=bam2df(ratab)
ratmockcount=bam2df(ratmock)
ratTEcount=bam2df(ratTE)

##total read numbers copied from command line wc -l reads.fq
readcounts=data.frame(sample=c('mock', 'infected', 'antibody'), numreads=c(1392516/4,1581172/4,1689296/4), sindbis=c(dim(mock)[1], dim(TE)[1], dim(ab)[1]), rat=c(dim(ratmockcount)[1], dim(ratTEcount)[1], dim(ratabcount)[1]))

pdf(file.path('~/Dropbox/Timplab_Data/sindbis/align_counts.pdf'), width=5, height=5)
print(grid.table(readcounts, rows=NULL))
dev.off()
