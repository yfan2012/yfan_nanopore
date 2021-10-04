library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/barnyard/barcode'
dbxdir='~/gdrive/mdr/barnyard'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)

prefix='210922_mdr_barnyard_speciescheck'

barcodefile=file.path(datadir, paste0(prefix, '_barcodes.txt'))
countsfile=file.path(datadir, paste0(prefix, '_barcodes_motifcounts.txt'))
fullbcinfo=read_tsv(barcodefile, col_name=bc_cols, na=c('None'))
fullbccounts=read_tsv(countsfile, col_name=bc_cols, na=c('None'))

nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
lowna=nacount[nacount<.25]
keepmotifs=names(lowna)

bcinfo=fullbcinfo %>%
    select(all_of(keepmotifs))


bccounts=fullbccounts %>%
    select(all_of(keepmotifs)) %>%
    filter(across(c(-readname, -chrname), ~.x>=3))
    
bcfilt=bcinfo %>%
    filter(readname %in% bccounts$readname)
    
ecoli=bcfilt %>%
    filter(chrname=='NC_000913.3') %>%
    filter(complete.cases(.))
staph=bcfilt %>%
    filter(chrname=='NC_007795.1') %>%
    filter(complete.cases(.))
strep=bcfilt %>%
    filter(chrname=='NZ_CP010450.1') %>%
    filter(complete.cases(.))
bacil=bcfilt %>%
    filter(chrname=='NC_000964.3') %>%
    filter(complete.cases(.))
    
counts=min(c(5000, dim(ecoli)[1], dim(staph)[1], dim(strep)[1], dim(bacil)[1]))

##all species
bcfiltsub=bind_rows(ecoli[1:counts,], staph[1:counts,], strep[1:counts,], bacil[1:counts,])
bcdata=bcfiltsub %>%
    select(-chrname, -readname)
bcumap=umap(bcdata)
plot=plot_umap(bcumap$layout, bcfiltsub$chrname, alpha=1, size=.4) +
    ggtitle('speciescheck')

##ecoli and staph
ecolistaph=bind_rows(ecoli[1:3000,], staph[1:3000,])
ecolistaphdata=ecolistaph %>%
    select(-chrname, -readname)
ecolistaphumap=umap(ecolistaphdata)
ecolistaphplot=plot_umap(ecolistaphumap$layout, ecolistaph$chrname, alpha=.8, size=.5) +
    ggtitle('ecoli and staph')

##ecoli bacil
ecolibacil=bind_rows(ecoli[1:counts,], bacil[1:counts,])
ecolibacildata=ecolibacil %>%
    select(-chrname, -readname)
ecolibacilumap=umap(ecolibacildata)
ecolibacilplot=plot_umap(ecolibacilumap$layout, ecolibacil$chrname, alpha=1, size=.4) +
    ggtitle('ecoli and bacil')

##ecoli strep
ecolistrep=bind_rows(ecoli[1:counts,], strep[1:counts,])
ecolistrepdata=ecolistrep %>%
    select(-chrname, -readname)
ecolistrepumap=umap(ecolistrepdata)
ecolistrepplot=plot_umap(ecolistrepumap$layout, ecolistrep$chrname, alpha=1, size=.4) +
    ggtitle('ecoli and strep')

##bacil strep
bacilstrep=bind_rows(bacil[1:counts,], strep[1:counts,])
bacilstrepdata=bacilstrep %>%
    select(-chrname, -readname)
bacilstrepumap=umap(bacilstrepdata)
bacilstrepplot=plot_umap(bacilstrepumap$layout, bacilstrep$chrname, alpha=1, size=.4) +
    ggtitle('bacil and strep')

##bacil staph
bacilstaph=bind_rows(bacil[1:counts,], staph[1:counts,])
bacilstaphdata=bacilstaph %>%
    select(-chrname, -readname)
bacilstaphumap=umap(bacilstaphdata)
bacilstaphplot=plot_umap(bacilstaphumap$layout, bacilstaph$chrname, alpha=1, size=.4) +
    ggtitle('bacil and staph')

##strep staph
strepstaph=bind_rows(strep[1:counts,], staph[1:counts,])
strepstaphdata=strepstaph %>%
    select(-chrname, -readname)
strepstaphumap=umap(strepstaphdata)
strepstaphplot=plot_umap(strepstaphumap$layout, strepstaph$chrname, alpha=1, size=.4) +
    ggtitle('strep and staph')



plotfile=file.path(dbxdir, '210922_mixes_speciescheck.pdf')
pdf(plotfile, w=13, h=7)
print(plot)
print(ecolistaphplot)
print(ecolibacilplot)
print(ecolistrepplot)
print(bacilstrepplot)
print(bacilstaphplot)
print(strepstaphplot)
dev.off()


