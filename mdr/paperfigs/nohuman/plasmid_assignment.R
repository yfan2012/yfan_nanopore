library(tidyverse)
library(RColorBrewer)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_nohuman'
datadir=file.path(projdir, 'paperfigs/nohuman')

dbxdir='~/gdrive/mdr/paperfigs/figs_nohuman'


##chrom to bin info, species info
chrombinsfile=file.path(datadir,'tigs2bins.tsv')
chrombins=read_tsv(chrombinsfile)

speciesfile=file.path(datadir, 'tigbin_species.tsv')
speciesinfo=read_tsv(speciesfile)


##grab plasmid info
tigplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native_nohuman.plasmidfinder.tsv')
plas_cols=c('file', 'seq', 'start', 'end', 'strand', 'gene', 'coverage', 'covmap', 'gaps', 'covfrac', 'ident', 'db', 'acc', 'prod', 'res')
tigplas=read_tsv(tigplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc)

##match to bin plasmid info
binplasfile=file.path(projdir, 'mdr/amr/200708_mdr_stool16native.hiC.plasmidfinder.tsv')
binplas=read_tsv(binplasfile, col_names=plas_cols, skip=1) %>%
    select(-file, -start, -end, -strand, -gaps, -coverage, -covmap, -covfrac, -db, -prod, -res, -acc)

allplas=full_join(tigplas, binplas, by='gene') %>%
    relocate(seq.x, .after=gene)
names(allplas)=c('gene', 'contig', 'tig.ident', 'bin', 'bin.ident')
plasinfo=rbind(allplas %>% na.omit, allplas %>% filter(is.na(contig)), allplas %>% filter(is.na(bin))) %>%
    filter(bin.ident>95 | is.na(bin.ident)) %>%
    filter(tig.ident>95 | is.na(tig.ident))
plasinfofile=file.path(datadir, 'plas_tigsbins.tsv')
write_tsv(plasinfo, plasinfofile)



##set up frequency info
methfreq=readRDS(file.path(datadir, 'methfreq.rds'))
invec=c('CAGAG','CCWGG', 'CMTCGAKG','CTCCAG', 'CTKVAG', 'GATC', 'GCGC', 'GCWGC', 'GGCC', 'GGNNCC', 'GGWCC', 'TCCGGA', 'GCCGGC', 'RGCGCY')
freqs=methfreq %>%
    rowwise() %>%
    filter(motif %in% invec) %>%
    spread(key=motif, value=freq)

plasfreqs=freqs %>%
    filter(chrom %in% tigplas$seq) %>%
    filter(complete.cases(.))

plasnear=NULL
for (i in tigplas$seq) {
    plasnearest=plascheck(freqs,i)
    plasnear=bind_rows(plasnear, plasnearest[1:10,])
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



####add kraken stuff
krakenfile=file.path(datadir, 'contig_id/kraken_standard/200708_mdr_stool16native_nohuman.out.short.txt')
kraken=read_tsv(krakenfile, col_names=c('class', 'tig', 'species'))
plasnear.kraken=plasnearbins %>%
    rowwise() %>%
    mutate(chroms2.kraken=kraken$species[kraken$tig==chroms2])

plasnearkrakencsv=file.path(dbxdir, 'plasmid_nearest_kraken.csv')
write_csv(plasnear.kraken, plasnearkrakencsv)
