library(tidyverse)
library(gplots)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

discodir='/mithril/Data/Nanopore/projects/methbin/disco'
datadir=file.path(discodir, 'barcode')
blastdir=file.path(discodir, 'blast_contigs')
dbxdir='~/gdrive/mdr/disco'


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


####identify contigs
contiginfofile=file.path(blastdir, 'classify_keys.csv')
tiglabs=read_csv(contiginfofile)


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



####heatmap with labels
tiginfolab=tiginfo %>%
    rowwise() %>%
    filter(chrname %in% tiglabs$tig) %>%
    mutate(label=tiglabs$name[which(tiglabs$tig==chrname)])

unscaledhm=as.matrix(tiginfolab[,2:8])
rownames(unscaledhm)=tiginfolab$label


distmat=dist(unscaledhm, method='euclidean')
hclustavg=hclust(distmat, method='average')

heatmappdf=file.path(dbxdir, 'contig_heatmap_unscaled.pdf')
pdf(heatmappdf, h=10, w=7)
hm=heatmap(unscaledhm, scale='none')
print(hm)
dev.off()

clustpdf=file.path(dbxdir, 'contig_hclust.pdf')
pdf(clustpdf, h=4, w=25)
plot(hclustavg)
dev.off()
