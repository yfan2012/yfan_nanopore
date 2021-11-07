library(tidyverse)

datadir='/mithril/Data/Nanopore/projects/methbin/zymo/contig_agg/20190809_zymo_control'
sumfile=file.path(datadir, '20190809_zymo_control.meth_report.txt')
dbxdir='~/gdrive/mdr/zymo'

sumfilecols=c('chr', 'pos', 'meth', 'unmeth')

summary=read_tsv(sumfile, col_names=sumfilecols) %>%
    filter(!grepl('tig', chr, fixed=TRUE)) %>%
    filter((meth+unmeth)>15) %>%
    mutate(methfrac=meth/(meth+unmeth))

fracpdf=file.path(dbxdir, 'megalodon_methfracs.pdf')
pdf(fracpdf, w=15, h=8)
plot=ggplot(summary, aes(x=methfrac, colour=chr, fill=chr, alpha=.3)) +
    geom_histogram() +
    scale_y_log10()+
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    facet_wrap(~chr) +
    theme_bw()
print(plot)
dev.off()
