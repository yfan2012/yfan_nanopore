library(ggplot2)
library(tidyverse)
library(R.utils)

##read in cov files, basic cov plot
samps=c('infected', 'mock', 'antibody')
##datadir='/scratch/groups/mschatz1/cpowgs/sindbis'
##samps=c('Antibody_2dpi', 'Antibody_3dpi', 'Sindbis_2dpi', 'Sindbis_3dpi')
datadir='/dilithium/Data/Nanopore/sindbis'
outdir='~/Dropbox/timplab_data/sindbis/cov/cov_plots'

allcov=NULL

##runstats=paste0('~/Dropbox/timplab_data/sindbis/runstats2.csv')
runstats=paste0('~/Dropbox/timplab_data/sindbis/runstats.csv')
yieldinfo=read_csv(runstats)

##annot_regions is a manually made regions file from the gb file in the dbox
gff=paste0(datadir, '/annot_regions.tsv')
regions=read_tsv(gff, col_names=c('prot', 'start', 'end')) %>%
    mutate(ymax=0) %>%
    mutate(colors=c('#EF9A9A', '#CE93D8', '#9FA8DA', '#81D4FA', '#80CBC4', '#C5E1A5', '#FFF59D', '#FFCC80', '#FFAB91'))

for (samp in samps) {
    covfile=paste0(datadir, '/',samp, '/cov/', samp, '.primary.cov') 
    fqfile=paste0(datadir, '/', samp, '/fqs/', samp, '.fq')
    yieldmb=yieldinfo$yield[yieldinfo$file==fqfile]/1000000
    cover=read_tsv(covfile, col_names=c('chr', 'pos', 'cov')) %>%
        mutate(sample=samp) %>%
        mutate(normcov=cov/yieldmb)
    sinvyieldmb=sum(cover$cov)/1000000
    cover=cover %>%
        mutate(sinvcov=cov/sinvyieldmb) %>%
        mutate(cov=cov/1000)
    
    allcov=bind_rows(allcov, cover)
    
    xwidth=max(cover$pos)-min(cover$pos)
    yheight=max(cover$cov)

    pdf(paste0(outdir,'/', samp, '.primary.cov.pdf'), width=20, height=5)
    barheight=-.1*max(cover$normcov)
    regions$ymin=barheight
    plot=ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
        geom_rect() +
        geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
        scale_fill_manual(values=regions$colors, labels=regions$prot) +
        geom_line(data=cover, inherit.aes=F, aes(x=pos, y=normcov))+
        xlab('Sindbis Genome Position') +
        ylab('Depth') +
        ggtitle(paste0(samp, ' Coverage (run yield normalized)')) +
        theme_bw()
    print(plot)
    barheight=-.1*max(cover$sinvcov)
    regions$ymin=barheight
    plot=ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
        geom_rect() +
        geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
        scale_fill_manual(values=regions$colors, labels=regions$prot) +
        geom_line(data=cover, inherit.aes=F, aes(x=pos, y=sinvcov))+
        xlab('Sindbis Genome Position') +
        ylab('Depth') +
        ggtitle(paste0(samp, ' Coverage (sindbis yield normalized)')) +
        theme_bw()
    print(plot)
    dev.off()
}

day1=rbind(allcov[allcov$sample=='antibody',], allcov[allcov$sample=='infected',])
day2=rbind(allcov[allcov$sample=='Antibody_2dpi',], allcov[allcov$sample=='Sindbis_2dpi',])
day3=rbind(allcov[allcov$sample=='Antibody_3dpi',], allcov[allcov$sample=='Sindbis_3dpi',])

pdf(paste0(outdir,'/all_primary_day1.normcov.pdf'), width=20, height=5)
##pdf(paste0(outdir,'/all_primary.normcov.pdf'), width=20, height=5)
barheight=-.1*max(cover$cov)
regions$ymin=barheight
print(ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
      geom_rect(show.legend=FALSE) +
      geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
      scale_fill_manual(values=regions$colors, labels=regions$prot) +
      geom_line(data=allcov, inherit.aes=F,  aes(x=pos, y=cov, colour=sample)) +
      xlab('Sindbis Genome Position') +
      ylab('Depth (kb)') +
      ggtitle(paste0('Coverage')) +
      theme_bw())
barheight=-.1*max(cover$normcov)
regions$ymin=barheight
print(ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
      geom_rect(show.legend=FALSE) +
      geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
      scale_fill_manual(values=regions$colors, labels=regions$prot) +
      geom_line(data=allcov, inherit.aes=F,  aes(x=pos, y=normcov, colour=sample)) +
      xlab('Sindbis Genome Position') +
      ylab('Depth (per mb yield)') +
      ggtitle(paste0('Coverage (run yield normalized)')) +
      theme_bw())
barheight=-.1*max(cover$sinvcov)
regions$ymin=barheight
print(ggplot(data=regions, mapping=aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
      geom_rect(show.legend=FALSE) +
      geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
      scale_fill_manual(values=regions$colors, labels=regions$prot) +
      geom_line(data=allcov, inherit.aes=F,  aes(x=pos, y=sinvcov, colour=sample)) +
      xlab('Sindbis Genome Position') +
      ylab('Depth (per mb sinv yield)') +
      ggtitle(paste0('Coverage (sindbis yield normalized)')) +
      theme_bw())

barheight=-.1*max(cover$normcov)
regions$ymin=barheight
print(ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
      geom_rect(show.legend=FALSE) +
      geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
      scale_fill_manual(values=regions$colors, labels=regions$prot) +
      geom_line(data=day2, inherit.aes=F,  aes(x=pos, y=normcov, colour=sample)) +
      xlab('Sindbis Genome Position') +
      ylab('Depth (per mb yield)') +
      theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1)) +
      ggtitle(paste0('Coverage (run yield normalized) Day 2')) +
      theme_bw() +
      theme(axis.text=element_text(size=24), axis.title=element_text(size=24)))
      

barheight=-.1*max(cover$normcov)
regions$ymin=barheight
print(ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
      geom_rect(show.legend=FALSE) +
      geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
      scale_fill_manual(values=regions$colors, labels=regions$prot) +
      geom_line(data=day3, inherit.aes=F,  aes(x=pos, y=normcov, colour=sample)) +
      xlab('Sindbis Genome Position') +
      ylab('Depth (per mb yield)') +
      ggtitle(paste0('Coverage (run yield normalized) Day 3')) +
      theme_bw()+
      theme(axis.text=element_text(size=24), axis.title=element_text(size=24)))

barheight=-.1*max(cover$normcov)
regions$ymin=barheight
print(ggplot(data=regions, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=prot), alpha=.3) +
      geom_rect(show.legend=FALSE) +
      geom_text(data=regions, aes(x=start+(end-start)/2, y=ymin+(ymax-ymin)/2, label=prot)) +
      scale_fill_manual(values=regions$colors, labels=regions$prot) +
      geom_line(data=day1, inherit.aes=F,  aes(x=pos, y=normcov, colour=sample)) +
      xlab('Sindbis Genome Position') +
      ylab('Depth (per mb yield)') +
      ggtitle(paste0('Coverage (run yield normalized) Day 1')) +
      theme_bw() +
      theme(axis.text=element_text(size=24), axis.title=element_text(size=24)))


dev.off()
