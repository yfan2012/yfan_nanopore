library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'
bcfile=file.path(datadir, '20190809_zymo_control_barcodes.txt')

##load data
bc_cols=c('readname', 'chrname', 'GATC', 'CCWGG', 'ATGCAT', 'GTCGAC', 'CTCCAG', 'CTKVAG')
bc=read_tsv(bcfile, col_names=bc_cols, na=c('', 'NA', 'None'))


##svd
bcdatasub=bc %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    group_by(chrname) %>%
    sample_n(20000) %>%
    ungroup()
bcdata=bcdatasub %>%
    select(-chrname, -readname)
bcdata[is.na(bcdata)]=0
#bcdata=t(bcdata)
s=svd(bcdata)
u=s$u
v=s$v
bcreduced4=u[,1:4]%*%v[1:4,]
bcreduced3=u[,1:3]%*%v[1:3,]

##scale all motifs to mean 0 and sd 1
bcscaled=apply(bcdata, 2, scale)
S=svd(bcscaled)
U=S$u
V=S$v
bcscaled4=U[,1:4]%*%V[1:4,]
bcscaled3=U[,1:3]%*%V[1:3,]

##umap
red4umap=umap(bcreduced4)
red3umap=umap(bcreduced3)

scale4umap=umap(bcscaled4)
scale3umap=umap(bcscaled3)


##plot
svdfile=file.path(dbxdir, 'umap_svd.pdf')
pdf(svdfile, h=7, w=13)
plot4=plot_umap(red4umap$layout, bcdatasub$chrname, c(-25,25), c(-25,25))
sep4=plot4 + facet_wrap(~label)
print(plot4)
print(sep4)
plot3=plot_umap(red3umap$layout, bcdatasub$chrname, c(-25,25), c(-25,25))
sep3=plot3 + facet_wrap(~label)
print(plot3)
print(sep3)
dev.off()

scalefile=file.path(dbxdir, 'umap_svd_scale.pdf')
pdf(scalefile, h=7, w=13)
plot4=plot_umap(scale4umap$layout, bcdatasub$chrname, c(-25,25), c(-25,25))
sep4=plot4 + facet_wrap(~label)
print(plot4)
print(sep4)
plot3=plot_umap(scale3umap$layout, bcdatasub$chrname, c(-25,25), c(-25,25))
sep3=plot3 + facet_wrap(~label)
print(plot3)
print(sep3)
dev.off()


