library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/qc/classify_plasmid_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo'
dbxdir='~/gdrive/mdr/zymo'
prefix='20190809_zymo_control_polished'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)
keyfile=file.path(datadir, 'medaka', 'consensus_key.csv')
key=read_csv(keyfile)

plasnames=c('Staphylococcus_aureus_plasmid1' ,'Escherichia_coli_plasmid')

##load meth calls from asm
barcodedatafile=file.path(datadir, 'barcode', prefix, '20190809_zymo_control_contigs_barcodes15.txt')
fullbcinfo=read_tsv(barcodedatafile, col_names=bc_cols, na=c('None'))
countsfile=file.path(datadir, 'barcode', prefix, '20190809_zymo_control_barcodes15_motifcounts.txt')
fullbccounts=read_tsv(countsfile, col_name=bc_cols, na=c('None'))
nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
lowna=nacount[nacount<.2]
keepmotifs=names(lowna)

##load meth calls from ref
refdatafile=file.path(datadir, 'barcode/20190809_zymo_control_barcodes15.txt')
fullrefinfo=read_tsv(refdatafile, col_names=bc_cols, na=c('None'))
refcountsfile=file.path(datadir, 'barcode/20190809_zymo_control_motifcounts15.txt')
refbccounts=read_tsv(refcountsfile, col_name=bc_cols, na=c('None'))
refna=colSums(is.na(fullrefinfo)/dim(fullrefinfo)[1])
refna=refna[refna<.2]
keepmotifs=names(refna)



####mason paper
##says CAGAG is a differential motif
refcounts=refbccounts %>%
    select(all_of(keepmotifs)) %>%
    filter(across(c(-readname, -chrname), ~ .x>=2))
reffilt=fullrefinfo %>%
    select(all_of(keepmotifs)) %>%
    filter(chrname %in% refcounts$chrname) %>%
    filter(complete.cases(.)) %>%
    filter(chrname %in% c('Escherichia_coli_chromosome', 'Salmonella_enterica_complete_genome')) %>%
    select(c(readname, chrname, CAGAG))


plotfile=file.path(dbxdir, 'ecoli_sal_explo.R')
ecolicolor=brewer.pal(3, 'Set2')[1]
salcolor=brewer.pal(3, 'Set2')[2]
pdf(plotfile, h=7, w=15)
CAGAGdensity=ggplot(reffilt, aes(x=CAGAG, colour=chrname, fill=chrname, alpha=.3)) +
    geom_histogram(data=subset(reffilt, chrname=='Escherichia_coli_chromosome'), fill=ecolicolor, colour=ecolicolor, alpha=.3) +
    geom_histogram(data=subset(reffilt, chrname=='Salmonella_enterica_complete_genome'), fill=salcolor, colour=salcolor, alpha=.3) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    xlim(-3,25) +
    theme_bw()
print(CAGAGdensity)
dev.off()

##investigate sal reads that have 0 CAGAG barcode score
salzero=reffilt %>%
    filter(chrname=='Salmonella_enterica_complete_genome') %>%
    filter(CAGAG==0)
salnonzero=reffilt %>%
    filter(chrname=='Salmonella_enterica_complete_genome') %>%
    filter(CAGAG!=0)

alignfile=file.path(datadir, 'align/20190809_zymo_control.paf')
pafcols=c('qname', 'qlen', 'qstart', 'qend', 'strand', 'rname', 'rlen', 'rstart', 'rend', 'matches', 'alen', 'mapq', 'tp', 'cm', 's1', 's2', 'dv', 'r1')
allpaf=read_tsv(alignfile, col_names=pafcols)
salzeropaf=allpaf %>%
    rowwise() %>%
    filter(qname %in% salzero$readname)
salnonzeropaf=allpaf %>%
    rowwise() %>%
    filter(qname %in% salnonzero$readname)

multialign=salzeropaf %>%
    filter(rname!='Salmonella_enterica_complete_genome')

salzerofiltpaf=salzeropaf %>%
    mutate(percmatch=matches/qlen) %>%
    filter(mapq==60)
salnonzerofiltpaf=salnonzeropaf %>%
    mutate(percmatch=matches/qlen) %>%
    filter(mapq==60)

##THE MOTIF DEFINING SALMONELLA FROM ECOLI IS NON PALINDROMIC WHICH IS WHY MY SHIT ISN'T DETECTING WELL. 



####see if the palindrome thing has been fixed
barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes20.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)
prefix='20190809_zymo_control'

plasnames=c('Staphylococcus_aureus_plasmid1' ,'Escherichia_coli_plasmid')

ref2datafile=file.path(datadir, 'barcode_v2', prefix, '20190809_zymo_control_barcodes20.txt')
fullref2info=read_tsv(ref2datafile, col_names=bc_cols, na=c('None'))
ref2countsfile=file.path(datadir, 'barcode_v2', prefix, '/20190809_zymo_control_motifcounts20.txt')
ref2bccounts=read_tsv(ref2countsfile, col_name=bc_cols, na=c('None'))
ref2na=colSums(is.na(fullref2info)/dim(fullref2info)[1])
ref2na=ref2na[ref2na<.2]
keepmotifs=names(ref2na)

ref2counts=ref2bccounts %>%
    select(all_of(keepmotifs)) %>%
    filter(across(c(-readname, -chrname), ~ .x>=2))
ref2filt=fullref2info %>%
    select(all_of(keepmotifs)) %>%
    filter(chrname %in% ref2counts$chrname) %>%
    filter(complete.cases(.)) %>%
    filter(chrname %in% c('Escherichia_coli_chromosome', 'Salmonella_enterica_complete_genome')) %>%
    select(c(readname, chrname, CAGAG))

plotfile=file.path(dbxdir, 'ecoli_sal_explo_v2.R')
ecolicolor=brewer.pal(3, 'Set2')[1]
salcolor=brewer.pal(3, 'Set2')[2]
pdf(plotfile, h=7, w=15)
CAGAGdensity=ggplot(ref2filt, aes(x=CAGAG, colour=chrname, fill=chrname, alpha=.3)) +
    geom_histogram(data=subset(ref2filt, chrname=='Escherichia_coli_chromosome'), fill=ecolicolor, colour=ecolicolor, alpha=.3) +
    geom_histogram(data=subset(ref2filt, chrname=='Salmonella_enterica_complete_genome'), fill=salcolor, colour=salcolor, alpha=.3) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    xlim(-3,25) +
    theme_bw()
print(CAGAGdensity)
dev.off()
