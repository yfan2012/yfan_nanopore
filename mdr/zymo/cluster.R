library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'


##load info
barcodefile=file.path(datadir, '20190809_zymo_control_barcodes.txt')
bc_cols=c('readname', 'chrname', 'GATC', 'CCWGG', 'ATGCAT','GTCGAC','CTCCAG','CTKVAG')
bcinfo=read_tsv(barcodefile, col_names=bc_cols, na=c('NA', 'None'))


##exclude yeast contigs which all start with 'tig' and plasmid
bc=bcinfo %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    group_by(chrname) %>%
    sample_n(10000) %>%
    ungroup()

    
umap_df <- function(info) {
    data=info %>%
        select(-readname, -chrname)
    data[is.na(data)]=0
    return(data)
}

elim_motifs  <- function(red, pdfname, x=NULL, y=NULL) {
    reddata=umap_df(red)
    redumap=umap(reddata)
    
    umappdf=file.path(dbxdir, pdfname)
    pdf(umappdf, h=7, w=13)
    if (missing(x)) {
        plot=plot_umap(redumap$layout, red$chrname)
    }else{
        plot=plot_umap(redumap$layout, red$chrname, x, y)
    }
    print(plot)
    sep=plot + facet_wrap(~label)
    print(sep)
    dev.off()
    
    return(redumap)
}


##count NAs
nacount=colSums(is.na(bcinfo)/dim(bcinfo)[1])


##try with 6 motifs
bcdata=umap_df(bc)
bcumap=umap(bcdata)

umappdf=file.path(dbxdir, 'umap_6motif.pdf')
pdf(umappdf, h=7, w=13)
plot=plot_umap(bcumap$layout, bc$chrname)
print(plot)
sep=plot + facet_wrap(~label)
print(sep)
plot=plot_umap(bcumap$layout, bc$chrname, c(-25,25), c(-25,25))
print(plot)
sep=plot + facet_wrap(~label)
print(sep)
dev.off()


##try with 5 motifs
red=bc %>%
    select(-GTCGAC)
elim_motif=elim_motifs(red, 'umap_5motif.pdf')
elim_motif_z=elim_motifs(red, 'umap_5motif_zoom.pdf', c(-25,25), c(-25,25))

##try with 4 motifs
red2=bc %>%
    select(-GTCGAC, -CTCCAG)
elim2_motif=elim_motifs(red2, 'umap_4motif.pdf')
elim2_motif_z=elim_motifs(red2, 'umap_4motif_zoom.pdf', c(-25,25), c(-25,25))

##try with 3 motifs
red3=bc %>%
    select(-GTCGAC, -CTCCAG, -ATGCAT)
elim3_motif=elim_motifs(red3, 'umap_3motif.pdf')
elim3_motif_z=elim_motifs(red3, 'umap_3motif_zoom.pdf', c(-25,25), c(-25,25))




##try throwing out reads
bcfilt=bcinfo %>%
    select(-GTCGAC, -CTCCAG)
bcnames=bcfilt[complete.cases(bcfilt),] %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    group_by(chrname) %>%
    sample_n(10000) %>%
    ungroup()
bcfilt=bcnames %>%
    select(-readname, -chrname)

filtumap=umap(bcfilt)

filtpdf=file.path(dbxdir, 'umap_readfilt.pdf')
pdf(filtpdf, h=7, w=13)
plot=plot_umap(filtumap$layout, bcnames$chrname)
sep=plot + facet_wrap(~label)
print(plot)
print(sep)
dev.off()
