library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='/mithril/Data/Nanopore/projects/methbin/disco/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/disco'


####contig based single read stuff
meta='MinION_JM3O_NAT_meta'

barcodesfile='~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt'
barcodes=read_tsv(barcodesfile, col_names=c('motif'))
bc_cols=c('readname', 'chrname', barcodes$motif)

bcfile=file.path(datadir, meta, paste0(meta, '_barcodes15.txt'))
bcinfo=read_tsv(bcfile, col_names=bc_cols, na=c('None'))

nas=colSums(is.na(bcinfo))/dim(bcinfo)[1]
keepnames=names(nas[nas<.3])


countsfile=file.path(datadir, meta, paste0(meta, '_motifcounts.txt'))
countsinfo=read_tsv(countsfile, col_names=bc_cols) %>%
    select(all_of(keepnames)) %>%
    filter(across(c(-readname, -chrname), ~ .x>3))
    
bcfilt=bcinfo %>%
    select(all_of(keepnames)) %>%
    filter(complete.cases(.)) %>%
    filter(readname %in% countsinfo$readname)

tigcounts=table(bcfilt$chrname)
tigsfilt=names(tigcounts[tigcounts>100])

tiginfo=bcfilt %>%
    select(-readname) %>%
    filter(chrname %in% tigsfilt) %>%
    group_by(chrname) %>%
    summarise(across(everything(), mean))


####try heatmap
tigscaled=scale(as.data.frame(tiginfo[,2:8]))
rownames(tigscaled)=tiginfo$chrname

heatmappdf=file.path(dbxdir, 'contig_heatmap.pdf')
pdf(heatmappdf, h=10, w=7)
hm=heatmap(tigscaled, scale='none')
print(hm)
dev.off()



####try umap
tigdata=tiginfo %>%
    select(-chrname)
contigumap=umap(tigdata)
umapdata=tibble(x=contigumap$layout[,1],
                y=contigumap$layout[,2],
                label=tiginfo$chrname)
contigumappdf=file.path(dbxdir, 'contig_umap.pdf')
pdf(contigumappdf, h=7, w=13)
plot=ggplot(umapdata, aes(x=x, y=y)) +
    geom_point(alpha=.2) +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
print(plot)
dev.off()


