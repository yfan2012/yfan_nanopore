library(tidyverse)
library(umap)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')

datadir='~/data/mdr/zymo/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'

barcodefile=file.path(datadir, '20190809_zymo_control_barcodes_filtered.txt')
bc_cols=c('readname', 'chrname', 'start', 'end', 'mapq', 'GATC', 'CCWGG', 'ATGCAT','GTCGAC','CTCCAG','CTKVAG')

bcinfo=read_tsv(barcodefile, col_names=bc_cols) %>%
    select(-start, -end, -mapq)

thresh=1.7
calls=call_part_avg(bcinfo, thresh)

orgs=unique(calls$chrname)
orgs=orgs[substring(orgs, 1,3)!='tig']

pops=calls %>%
    group_by(chrname) %>%
    summarise(pops=as.numeric(table(barcode)), barcode=names(table(barcode)))

orgbcpdf=file.path(dbxdir, 'zymo_org_bc.pdf')
pdf(orgbcpdf, w=19, h=9)
for (i in orgs) {
    orgpops=pops %>%
        filter(chrname==i) %>%
        rename(samp=chrname)
    plot=plot_pops(orgpops, i)
    print(plot)
}
dev.off()



####try umap
bacbcinfo=bcinfo %>%
    rowwise() %>%
    filter(chrname %in% orgs)

bcdata=bacbcinfo %>%
    select(-readname, -chrname)

bcumap=umap(bcdata)

plot_umap <- function(layout, labels) {
    data=tibble(x=layout[,1],
                y=layout[,2],
                label=labels)
    plot=ggplot(data, aes(x=x, y=y, colour=label, alpha=.1)) +
        geom_point() +
        scale_colour_brewer(palette='Set2') +
        theme_bw()
    return(plot)
}

umappdf=file.path(dbxdir, 'umap.pdf')
pdf(umappdf, h=7, w=13)
print(plot_umap(bcumap$layout, bacbcinfo$chrname))
dev.off()
