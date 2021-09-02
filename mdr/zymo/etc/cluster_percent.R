library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')


datadir='/mithril/Data/Nanopore/projects/methbin/zymo/contig'
dbxdir='~/gdrive/mdr/zymo'
prefix='20190809_zymo_control'

methfile=file.path(datadir, paste0(prefix, '_contig15.txt'))
covfile=file.path(datadir, paste0(prefix, '_contig15_cov.txt'))

meth=read_tsv(methfile) %>%
    filter(!grepl('tig', chrname, fixed=TRUE))

cov=read_tsv(covfile) %>%
    filter(!grepl('tig', chrname, fixed=TRUE))

bcinfo=as_tibble(meth[,-1]/cov[,-1]) 

nacount=colSums(is.na(bcinfo))/dim(bcinfo)[1]
keepnames=names(nacount[nacount<.2])

bcfilt=bcinfo %>%
    select(all_of(keepnames)) %>%
    mutate(chrname=meth$chrname) %>%
    filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>%
    filter(complete.cases(.))

tiginfo=bcfilt %>%
    group_by(chrname) %>%
    summarise(across(everything(), mean)) %>%
    ungroup()
clusttig=as.matrix(tiginfo %>% select(-chrname))[1:9,]
rownames(clusttig)=tiginfo$chrname[1:9]
tigdata=as.data.frame(tiginfo[1:9,-1])
scaletig=scale(tigdata)
rownames(scaletig)=tiginfo$chrname[1:9]


##get distance plot
##http://www.opiniomics.org/you-probably-dont-understand-heatmaps/
heatmappdf=file.path(dbxdir, 'heatmap_contig_calls_filt.pdf')
pdf(heatmappdf, h=7, w=13)
hm=heatmap(scaletig, scale='none')
noscale=heatmap(clusttig, scale='none')
tiginfodist=plot(hclust(dist(clusttig)))
distplot=plot(hclust(dist(scaletig)))
print(hm)
print(noscale)
print(distplot)
print(tiginfodist)
dev.off()


##look at why the staph plasmid clusters well with umap but doesn't with hierarchical
plascov=cov %>%
    filter(grepl('plasmid', chrname, fixed=TRUE))

staphcov=plascov %>%
    filter(grepl('Staphylococcus_aureus_plasmid1', chrname, fixed=TRUE))
staphchrcov=cov %>%
    filter(grepl('Staphylococcus_aureus_chromosome', chrname, fixed=TRUE))


##try umap with this stuff
bcchr=bcfilt %>%
    filter(!grepl('tig', chrname, fixed=TRUE)) %>%
    filter(!grepl('plasmid', chrname, fixed=TRUE)) %>%
    group_by(chrname) %>%
    sample_n(10000)
bcplas=bcfilt %>%
    filter(grepl('plasmid', chrname, fixed=TRUE))
allfilt=bind_rows(bcplas, bcchr)

alldata=allfilt %>%
    select(-chrname)

allumap=umap(alldata)
embed=tibble(x=allumap$layout[,1],
             y=allumap$layout[,2],
             label=allfilt$chrname) %>%
    mutate(colour=case_when(!grepl('plasmid', label, fixed=TRUE) ~ label, TRUE ~ 'plasmid')) %>%
    mutate(shape=case_when(grepl('plasmid', label, fixed=TRUE) ~ label, TRUE ~ 'chr'))

embedplas=embed %>%
    filter(colour=='plasmid') %>%
    select(-label)
embedchr=embed %>%
    filter(shape=='chr')


mycolors=c(brewer.pal(8, 'Set2'), '#000000')
myshapes=c(3,18, 1, 2)
plasfile=file.path(dbxdir, 'plasmid_shown_percent.pdf')
pdf(plasfile, h=9, w=13)
plot=ggplot(embedchr, aes(x=x, y=y, colour=colour)) +
    geom_point(alpha=.2, size=.1) +
    scale_colour_manual(values=mycolors) +
    ylim(-10,10) +
    xlim(-10,10) +
    theme_bw()
mainplot=plot +
    geom_point(data=embedplas, aes(x=x, y=y, shape=shape), inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes) +
    ylim(-10,10) +
    xlim(-10,10) +
    theme(legend.position = 'none')
print(mainplot)
sep=plot +
    facet_wrap(~label) +
    geom_point(data=embedplas, aes(x=x, y=y, shape=shape), size=.5, alpha=.4, inherit.aes=FALSE) +
    scale_shape_manual(values=myshapes) +
    ylim(-10,10) +
    xlim(-10,10) +
    theme(legend.position = "none")
print(sep)
dev.off()
