library(DESeq2)
library(tximport)
library(readr)
library(tidyverse)
library(cowplot)

datadir='/pym/Data/Nanopore/projects/prolificans/rna'
dbxdir='~/Dropbox/timplab_data/prolificans/rna'

conds=c('AMB', 'ITR', 'ND', 'VRC')
reps=c('1', '2', '3')
samps=tibble(reps=rep(reps, 4), cond=rep(conds, each=3)) %>%
    mutate(run=paste0(cond, reps))

##from vignette
files=file.path(datadir, 'salmon', paste0('5317_',samps$run), 'quant.sf')
names(files)=samps$cond

transmapfile=file.path(datadir, 'trinity', 'Trinity.fasta.gene_trans_map')
tx2gene=read_tsv(transmapfile, col_names=c('GENEID', 'TXNAME')) %>%
    relocate(TXNAME, GENEID)

txi=tximport(files, type="salmon", tx2gene=tx2gene)
dds=DESeqDataSetFromTximport(txi, colData=samps, design= ~ cond)
dds=DESeq(dds)
res=results(dds)

amb=data.frame(results(dds, contrast=c("cond","ND","AMB")))
itr=data.frame(results(dds, contrast=c("cond","ND","ITR")))
vrc=data.frame(results(dds, contrast=c("cond","ND","VRC")))


volcano_plot <- function(df, title) {
    ##volcano plot from deseq results
    tib=as_tibble(df) %>%
        mutate(gene=rownames(df)) %>%
        drop_na() %>%
        mutate(logpadj=-log10(padj))
    plot=ggplot(tib, aes(x=log2FoldChange, y=logpadj, alpha=.05)) +
        geom_point() +
        ggtitle(title) +
        xlab('log2 Fold Change') +
        xlim(-12,12) +
        ylab('-log10 adjusted p value') +
        ylim(0, 175) +
        theme_bw()
    return(plot)
}

ambplot=volcano_plot(amb, 'ND vs AMB')
itrplot=volcano_plot(itr, 'ND vs ITR')
vrcplot=volcano_plot(vrc, 'ND vs VRC')

volcanopdf=file.path(dbxdir, 'volcano.pdf')
pdf(volcanopdf, w=20, h=6)
plot_grid(ambplot, itrplot, vrcplot, ncol=3, align='h')
dev.off()


goinfotsv=file.path(datadir, 'transcriptome', 'st5317.goterms.tsv')
goinfo=read_tsv(goinfotsv, col_names=c('gene', 'goterms'))
desinfotsv=file.path(datadir, 'transcriptome', 'st5317.description.tsv')
desinfo=read_tsv(desinfotsv, col_names=c('gene', 'description'))

topgenes <- function(df, samp, goinfo, desinfo) {
    ##get most significant and highest change genes
    tib=as_tibble(df) %>%
        mutate(gene=rownames(df)) %>%
        drop_na() %>%
        mutate(logpadj=-log10(padj)) %>%
        mutate(importance=logpadj*abs(log2FoldChange)) %>%
        arrange(-importance) %>%
        mutate(samp=samp) %>%
        rowwise() %>%
        mutate(des=desinfo$description[desinfo$gene==gene])

    ##not sure why mutate/case_when didn't work in my testing
    go=c()
    for (i in 1:length(tib$gene)) {
        if (tib$gene[i] %in% goinfo$gene) {
            go[i]=goinfo$goterms[goinfo$gene==tib$gene[i]]
        }else{
            go[i]='unknown'
        }
    }
    tib$go=go
            
    top=tib[1:30,]
    return(top)
}

ambtop=topgenes(amb, 'AMB', goinfo, desinfo)
itrtop=topgenes(itr, 'ITR', goinfo, desinfo)
vrctop=topgenes(vrc, 'VRC', goinfo, desinfo)

topgenes=bind_rows(ambtop, itrtop, vrctop)
topcsv=file.path(dbxdir, 'top_diff_exp.csv')
write_csv(topgenes, topcsv)


