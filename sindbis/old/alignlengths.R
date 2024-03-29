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


##plot sindbis primary alignment lengths
abbam=file.path('/dilithium/Data/Nanopore/sindbis/antibody/align/antibody.primary.sorted.bam')
mockbam=file.path('/dilithium/Data/Nanopore/sindbis/mock/align/mock.primary.sorted.bam')
TEbam=file.path('/dilithium/Data/Nanopore/sindbis/infected/align/infected.primary.sorted.bam')
ab=bam2df(abbam)
ab$samp='antibody'
mock=bam2df(mockbam)
mock$samp='mock'
TE=bam2df(TEbam)
TE$samp='infected'
allsamps=rbind(ab, mock, TE)
v3samps=rbind(ab, TE)

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

pdf(file.path('~/Dropbox/timplab_data/sindbis/align/sindbis_pilot_align_lengths.pdf'), width=11, height=8.5)
print(ggplot(v3samps, aes(x=samp, y=rlength, fill=samp)) +
      geom_violin(scale="width") +
      ggtitle('Sindbis Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw() +
      theme(axis.text=element_text(size=24), axis.title=element_text(size=24)))
dev.off()


##look at the spliced alignments to the genome to count the rat rna
ratab=file.path('/dilithium/Data/Nanopore/sindbis/antibody/align/antibody.rat.splicealn.primary.sorted.bam')
ratmock=file.path('/dilithium/Data/Nanopore/sindbis/mock/align/mock.rat.splicealn.primary.sorted.bam')
ratTE=file.path('/dilithium/Data/Nanopore/sindbis/infected/align/infected.rat.splicealn.primary.sorted.bam')

ratabcount=bam2df(ratab)
ratmockcount=bam2df(ratmock)
ratTEcount=bam2df(ratTE)


##make a table of the counts of sindbis rna vs rat rna
##total read numbers copied from command line wc -l reads.fq
readcounts=data.frame(sample=c('mock', 'infected', 'antibody'), numreads=c(1392516/4,1581172/4,1689296/4), sindbis=c(dim(mock)[1], dim(TE)[1], dim(ab)[1]), rat=c(dim(ratmockcount)[1], dim(ratTEcount)[1], dim(ratabcount)[1]))

pdf(file.path('~/Dropbox/Timplab_Data/sindbis/align_counts.pdf'), width=5, height=5)
print(grid.table(readcounts, rows=NULL))
dev.off()



##plot alignment lengths of transcriptome
scriptab=file.path('/dilithium/Data/Nanopore/sindbis/antibody/align/antibody.rat.transcriptaln.primary.sorted.bam')
scriptmock=file.path('/dilithium/Data/Nanopore/sindbis/mock/align/mock.rat.transcriptaln.primary.sorted.bam')
scriptTE=file.path('/dilithium/Data/Nanopore/sindbis/infected/align/infected.rat.transcriptaln.primary.sorted.bam')
scriptabcount=bam2df(scriptab)
scriptabcount$samp='antibody'
scriptmockcount=bam2df(scriptmock)
scriptmockcount$samp='mock'
scriptTEcount=bam2df(scriptTE)
scriptTEcount$samp='infected'

scriptallsamps=rbind(scriptabcount, scriptmockcount, scriptTEcount)


pdf(file.path('~/Dropbox/Timplab_Data/sindbis/rat_align_lengths.pdf'), width=11, height=8.5)
print(ggplot(scriptallsamps, aes(x=samp, y=rlength, fill=samp)) +
      geom_violin(scale="width") +
      ggtitle('Rat Transcriptome Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
print(ggplot(scriptallsamps, aes(x=samp, y=rlength, fill=samp)) +
      geom_boxplot() +
      ggtitle('Rat Transcriptome Primary Alignment Lengths') +
      xlab('Sample') +
      ylab('Read Length') +
      theme_bw())
dev.off()
