library(tidyverse)
library(RColorBrewer)

dbxdir='~/gdrive/mdr/zymo'
datadir='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite/bismark'
labels=c('bsubtilis', 'ecoli', 'efaecalis', 'lmonocytogenes', 'nasa_Ecoli_K12', 'paeruginosa', 'saureus', 'senterica')
cxcols=c('chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'seq')


####sanity check: plot %methylation distributions
pdffile=file.path(dbxdir, 'bisulfite_methfracs.pdf')
pdf(pdffile, h=9, w=15)
mycolors=colorRampPalette(brewer.pal(8, 'Set2'))(11)
for (i in labels) {
    cxfile=file.path(datadir, i, paste0(i, '_1_bismark_bt2_pe.CX_report.txt'))
    cx=read_tsv(cxfile, col_names=cxcols) %>%
        filter(!grepl('tig', chr, fixed=TRUE)) %>%
        filter((meth+unmeth)>15) %>%
        mutate(methfrac=meth/(meth+unmeth))
    plot=ggplot(cx, aes(x=methfrac, colour=chr, fill=chr, alpha=.3)) +
        geom_histogram() +
        ggtitle(i) +
        scale_y_log10()+
        scale_colour_manual(values=mycolors) +
        scale_fill_manual(values=mycolors) +
        facet_wrap(~chr) +
        theme_bw()
    print(plot)
}
dev.off()

    
