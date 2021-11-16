library(tidyverse)
library(RColorBrewer)

prefix='20190809_zymo_control'
projdir='/mithril/Data/Nanopore/projects/methbin/zymo'
datadir=file.path(projdir, 'contig_agg', prefix)
sumfile=file.path(datadir, paste0(prefix, '.meth_report.txt'))
dbxdir='~/gdrive/mdr/zymo'


####plot methfracs from aggregation of meth text file
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


####write out summary of putative motifs
barcodefile='~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt'
bc=read_tsv(barcodefile, col_names=c('motif'))
bc_cols=c('chr', bc$motif)

motifinfofile=file.path(datadir, '20190809_zymo_control_barcodes50.csv')
motifinfo=read_csv(motifinfofile, col_names=bc_cols, na=c("", "None"))

summary=c()
for (i in 1:dim(motifinfo)[1]) {
    row=motifinfo[i,][-1]
    row[is.na(row)]=0
    a=row[as.vector(row>.8)]
    if (length(a)>0) {
        motifs=names(a)
        str=''
        for (j in 1:length(a)) {
            str=paste0(str, '\t', motifs[j], ',', as.character(a[j]))
        }
        info=paste0(motifinfo[i,1],str)
    }else{
        info=paste0(motifinfo[i,1])
    }
    summary=c(summary, info)
}

summaryfile=file.path(datadir, '20190809_zymo_control_motifs.txt')
write(summary, summaryfile)


####plot meth fracs from megalodon aggregation
megacols=c('chrom', 'start', 'end', 'name', 'score', 'strand', 'tstart', 'tend', 'color', 'cov', 'methfrac')
mycolors=colorRampPalette(brewer.pal(8, 'Set2'))(9)

cmethfile=file.path(projdir, 'megalodon', prefix, 'modified_bases.5mC.bed')
cmeth=read_tsv(cmethfile, col_names=megacols) %>%
    filter(!grepl('tig', chrom, fixed=TRUE)) %>%
    filter(cov>15)

amethfile=file.path(projdir, 'megalodon', prefix, 'modified_bases.6mA.bed')
ameth=read_tsv(amethfile, col_names=megacols) %>%
    filter(!grepl('tig', chrom, fixed=TRUE)) %>%
    filter(cov>15)

allmeth=bind_rows(ameth, cmeth)

aggmegapdf=file.path(dbxdir, 'agg_by_mega_methfracs.pdf')
pdf(aggmegapdf, w=15, h=8)
plot=ggplot(allmeth, aes(x=methfrac, colour=chrom, fill=chrom, alpha=.3)) +
    geom_histogram() +
    scale_y_log10()+
    ggtitle('All meth') +
    scale_colour_manual(values=mycolors) +
    scale_fill_manual(values=mycolors) +
    facet_wrap(~chrom) +
    theme_bw()
print(plot)
aplot=ggplot(ameth, aes(x=methfrac, colour=chrom, fill=chrom, alpha=.3)) +
    geom_histogram() +
    scale_y_log10()+
    ggtitle('A meth') +
    scale_colour_manual(values=mycolors) +
    scale_fill_manual(values=mycolors) +
    facet_wrap(~chrom) +
    theme_bw()
print(aplot)
cplot=ggplot(cmeth, aes(x=methfrac, colour=chrom, fill=chrom, alpha=.3)) +
    geom_histogram() +
    scale_y_log10()+
    ggtitle('C meth') +
    scale_colour_manual(values=mycolors) +
    scale_fill_manual(values=mycolors) +
    facet_wrap(~chrom) +
    theme_bw()
print(cplot)
dev.off()
