library(tidyverse)
library(RColorBrewer)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/figs'


##chrom to bin info
chrombinsfile='/mithril/Data/Nanopore/projects/methbin/paperfigs/contig_level/tigs2bins.tsv'
chrombins=read_tsv(chrombinsfile)


##grab plasmid info
tigplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.plasmidfinder.tsv')
plas_cols=c('file', 'seq', 'start', 'end', 'strand', 'gene', 'coverage', 'covmap', 'gaps', 'covfrac', 'ident', 'db', 'acc', 'prod', 'res')
tigplas=read_tsv(tigplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc)


##set up frequency info
exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))
##methfreq=readRDS(file.path(datadir, 'clin_methfreq3.rds'))
freqs=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec) %>%
    spread(key=motif, value=freq)

plasfreqs=freqs %>%
    filter(chrom %in% tigplas$seq) %>%
    filter(complete.cases(.))

plasnear=NULL
for (i in tigplas$seq) {
    plasnearest=plascheck(freqs,i)
    plasnear=bind_rows(plasnear, plasnearest[1:20,])
}
plasnear=unique(plasnear[complete.cases(plasnear),]) %>%
    select(-rounded)

plasnearbins=NULL
for (i in 1:dim(plasnear)[1]) {
    info=plasnear[i,]
    if (info$chroms %in% chrombins$rname) {
        info$bin1=chrombins$bin[chrombins$rname==info$chroms]
    } else {
        info$bin1='unknown'
    }

    if (info$chroms2 %in% chrombins$rname) {
        info$bin2=chrombins$bin[chrombins$rname==info$chroms2]
    } else {
        info$bin2='unknown'
    }
    plasnearbins=bind_rows(plasnearbins, info)
}

plasnearbinscsv=file.path(dbxdir, 'plasmid_nearest.csv')
write_csv(plasnearbins, plasnearbinscsv)

nearest=plasnearbins %>%
    filter(bin1!='unknown') %>%
    filter(bin2!='unknown') %>%
    rowwise() %>%
    filter(chroms %in% chrombins$rname) %>%
    group_by(chroms) %>%
    filter(row_number()==1) %>%
    arrange(-nummoitfs)


plascounts=plasnearbins %>%
    filter(bin1!='unknown') %>%
    group_by(chroms) %>%
    summarise(unk=sum(bin2=='unknown'), cor=sum(bin1==bin2), incor=sum(bin1!=bin2)-sum(bin2=='unknown')) %>%
    rowwise() %>%
    mutate(nummotifs=nearest$nummoitfs[nearest$chroms==chroms]) %>%
    gather("status", "counts", -c("nummotifs", "chroms")) %>%
    mutate(chromlabel=paste0(chroms, ': ', as.character(nummotifs), ' motifs')) %>%
    filter(status!='unk') %>%
    arrange(-nummotifs)

plascounts=within(plascounts, chromlabel <- factor(chromlabel, levels=unique(plascounts$chromlabel)))
colors=rev(brewer.pal(8, 'Set2')[4:5])
names(colors)=c('cor', 'incor')

plasdists=plasnearbins %>%
    filter(bin1!='unknown') %>%
    mutate(status=case_when(bin1==bin2 ~ 'cor',
                            bin2=='unknown' ~ 'unk',
                            TRUE ~ 'incor')) %>%
    rowwise() %>%
    mutate(nummotifs=nearest$nummoitfs[nearest$chroms==chroms]) %>%
    select(-c(chroms2, bin1, bin2)) %>%
    mutate(chromlabel=paste0(chroms, ': ', as.character(nummotifs), ' motifs'))

plasdists=within(plasdists, chromlabel <- factor(chromlabel, levels=unique(plascounts$chromlabel)))
distcolors=c(rev(brewer.pal(8, 'Set2')[4:5]), '#B0A1A0')
names(distcolors)=c('cor', 'incor', 'unk')


plasplotpdf=file.path(dbxdir, 'plasmid_assignments.pdf')
pdf(plasplotpdf, h=6, w=11)
plasplot=ggplot(plascounts, aes(x=status, y=counts, colour=status, fill=status, alpha=.4)) +
    geom_bar(stat='identity') +
    facet_wrap(~chromlabel, ncol=6) +
    ylim(0,20) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    theme_bw()
print(plasplot)
##https://stackoverflow.com/questions/49330742/change-y-axis-of-dot-plot-to-reflect-actual-count-using-geom-dotplot
##for adjusting y axis height
distplot=ggplot(plasdists, aes(x=dist, colour=status, fill=status, alpha=.8)) +
    geom_dotplot(stackgroups=T, method='histodot', stackratio=1) +
    facet_wrap(~chromlabel, ncol=5) +
    coord_fixed(ratio=(1.5/30)*18) +
    ##scale_y_continuous(NULL, breaks = NULL) +
    scale_color_manual(values=distcolors) +
    scale_fill_manual(values=distcolors) +
    theme_bw()
plot(distplot)
dev.off()



    
 





