library(tidyverse)
library(RColorBrewer)
library(Biostrings)

dbxdir='~/gdrive/mdr/zymo'
projdir='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite'
datadir=file.path(projdir, 'bismark')
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

    
####report motifs
##barcodefile='~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt'
barcodefile='~/Code/yfan_nanopore/mdr/rebase/barcodes51.txt'
bc=read_tsv(barcodefile, col_names=c('motif'))
bc_cols=c('chr', bc$motif)

motifinfofile=file.path(projdir, 'motifcalls', '51_barcodes.csv')
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

##summaryfile=file.path(projdir, 'motifcalls', 'meth_summary.txt')
summaryfile=file.path(projdir, 'motifcalls', 'meth_summary_51_barcodes.txt')
write(summary, summaryfile)



