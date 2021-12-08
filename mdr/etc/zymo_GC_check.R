library(tidyverse)
library(multidplyr)

cluster=new_cluster(12)
cluster_library(cluster, 'tidyverse')

dbxdir='~/gdrive/mdr/etc'
projdir='/mithril/Data/Nanopore/projects/methbin'
datadir=file.path(projdir, 'etc/CG_model')

methcall_cols=c('chrom', 'pos', 'strand', 'prob', 'motif', 'base', 'meth')
methcallsfile=file.path(datadir, 'curate_calls.csv')


##read motif info
cmethfile=file.path(projdir, 'zymo/truth/bisulfite/zymo_cmeth.csv')
cmeth=read_csv(cmethfile)
amethfile=file.path(projdir, 'zymo/truth/pacbio/zymo_ameth.csv')
ameth=read_csv(amethfile)
methinfo=bind_rows(cmeth,ameth) %>%
    group_by(motif, pos) %>%
    summarise(num=n())

nametochrfile='~/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
chrlabels=read_table2(nametochrfile, col_names=c('chr', 'label'))

meth=read_csv(methcallsfile, col_names=methcall_cols) %>%
    filter(!grepl('tig0', chrom, fixed=TRUE)) %>%
    group_by(chrom, pos, strand, motif) %>%
    summarise(methnum=sum(meth=='m'), umethnum=sum(meth=='u')) %>%
    mutate(methfrac=methnum/(methnum+umethnum))

findMethFreq <- function(x) {
    ##collapses called meth into methfreq.
    ##x is meth df above, grouped by chrom and motif

    motiflen=nchar(x$motif[1])

    motifpos=NULL
    for (i in x$pos) {
        diffs=abs(x$pos-i)
        motifgroup=x[diffs<=motiflen,] %>%
            arrange(-methfrac) %>%
            mutate(totcalls=methnum+umethnum) %>%
            arrange(-totcalls)
        motifpos=bind_rows(motifpos, motifgroup[1,])
    }

    motifpos=unique(motifpos) %>%
        select(-totcalls)
    return(motifpos)
}
cluster_copy(cluster, 'findMethFreq')

methgroups=meth %>%
    group_by(chrom, motif) %>%
    partition(cluster)
methfreqs=methgroups %>%
    do(findMethFreq(.)) %>%
    collect() %>%
    ungroup()
cpgs=methfreqs %>%
    filter(motif=='GCCGGC')


cgpdf=file.path(dbxdir, 'cg_model_megalodon.pdf')
pdf(cgpdf, h=9, w=13)
plot=ggplot(cpgs, aes(x=chrom, y=methfrac, colour=chrom, fill=chrom, alpha=.2)) +
    geom_boxplot() +
    ggtitle('GCCGGC') +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(plot)
dev.off()



labels=c('bsubtilis', 'ecoli', 'efaecalis', 'lmonocytogenes', 'paeruginosa', 'saureus', 'senterica')
cxcols=c('chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'seq')

allcx=NULL
for (label in labels) {
    cxfile=file.path(projdir, 'zymo/truth/bisulfite/bismark', label, paste0(label, '_1_bismark_bt2_pe.CX_report.txt'))
    cx=read_tsv(cxfile, col_names=cxcols) %>%
        mutate(total=meth+unmeth)
    allcx=bind_rows(allcx, cx)
}

methcx=allcx %>%
    rowwise() %>%
    filter(!grepl('tig', chr, fixed=TRUE)) %>%    
    group_by(chr, pos) %>%
    top_n(1, total)
