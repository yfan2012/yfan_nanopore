library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'

##load info
bc_cols=c('readname', 'chrname', 'GATC', 'CCWGG', 'ATGCAT','GTCGAC','CTCCAG','CTKVAG')
countsfile=file.path(datadir, '20190809_zymo_control_motifcounts.txt')
bccounts=read_tsv(countsfile, col_names=bc_cols, na=c('NA', 'None'))

barcodefile=file.path(datadir, '20190809_zymo_control_barcodes.txt')
bcinfo=read_tsv(barcodefile, col_names=bc_cols, na=c('NA', 'None'))

##count NA
nacount=colSums(is.na(bcinfo)/dim(bcinfo)[1])

##select for low na motifs
countsfiltered=tibble(mins=seq(1,10,1),
                      counts=rep(0, 10),
                      filtered=rep(0,10))
for (i in 1:10) {
    filteredcounts=bccounts %>%
        filter(across(c(-readname, -chrname, -GTCGAC, -CTCCAG), ~ .x>=i))
    countsfiltered$filtered[i]=dim(filteredcounts)[1]
    counts=bccounts %>%
        filter(across(c(-readname, -chrname), ~ .x>=i))
    countsfiltered$counts[i]=dim(counts)[1]
}
plotcounts=gather(countsfiltered, 'cond', 'counts', -mins)

countpdf=file.path(dbxdir, 'readcounts.pdf')
pdf(countpdf, h=8, w=11)
ggplot(plotcounts, aes(x=mins, y=counts, colour=cond, fill=cond)) +
    geom_line() +
    geom_point(size=3) +
    ggtitle('Filtered read counts') +
    xlab('Min number of motifs') +
    ylab('number of reads surviving') +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
dev.off()


##filter for more than 3 motifs
countsfilt=bccounts %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    select(-GTCGAC, -CTCCAG) %>%
    filter(across(c(-readname, -chrname), ~ .x>=3)) %>%
    group_by(chrname) %>%
    sample_n(10000) 
filtinfo=bcinfo %>%
    rowwise() %>%
    filter(readname %in% countsfilt$readname)

filtdata=countsfilt %>%
    ungroup() %>%
    select(-readname, -chrname)
filtumap=umap(filtdata)

umapfiltfile=file.path(dbxdir, 'filt_umap.pdf')
pdf(umapfiltfile, h=7, w=11)
plot=plot_umap(filtumap$layout, countsfilt$chrname)
print(plot)
sep=plot + facet_wrap(~label)
print(sep)
dev.off()

