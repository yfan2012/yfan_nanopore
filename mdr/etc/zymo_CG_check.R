library(tidyverse)
library(multidplyr)
library(Biostrings)

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
    filter(motif=='GCCGGC') %>%
    rename(chrom='chr')


cgpdf=file.path(dbxdir, 'cg_model_megalodon.pdf')
pdf(cgpdf, h=9, w=13)
plot=ggplot(cpgs, aes(x=chrom, y=methfrac, colour=chrom, fill=chrom, alpha=.2)) +
    geom_boxplot() +
    ggtitle('GCCGGC') +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(plot)
dev.off()




####bisulf stuff
keyfile='~/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
keycols=c('chrom', 'label')
key=read_table2(keyfile, col_names=keycols)

labels=c('bsubtilis', 'ecoli', 'efaecalis', 'lmonocytogenes', 'paeruginosa', 'saureus', 'senterica')
cxcols=c('chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'seq')
allcx=NULL
for (label in labels) {
    cxfile=file.path(projdir, 'zymo/truth/bisulfite/bismark', label, paste0(label, '_1_bismark_bt2_pe.CX_report.txt'))
    chrlist=key$chrom[key$label==label]
    cx=read_tsv(cxfile, col_names=cxcols) %>%
        rowwise() %>%
        filter(chr %in% chrlist) %>%
        mutate(total=meth+unmeth)
    allcx=bind_rows(allcx, cx)
}

methbisulf=allcx %>%
    mutate(bisulf=meth/total) 




####filter out appropriate 
reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
ref=readDNAStringSet(reffile, format='fasta')
chrs=names(ref)[!grepl('tig', names(ref), fixed=TRUE)]
methcompare=NULL
methcomparefixed=NULL
for (i in chrs) {
    print(i)
    chromname=strsplit(i, split=' ', fixed=TRUE)[[1]][1]
    starts=start(vmatchPattern('GCCGGC', ref[i])[[i]])
    bisulfpos=c(starts+2, starts+3)
    bisulf=methbisulf %>%
        filter(chr==chromname) %>%
        rowwise() %>%
        filter(pos %in% bisulfpos) %>%
        mutate(motif='GCCGGC')
    mega=meth %>%
        filter(chrom==chromname) %>%
        filter(motif=='GCCGGC')
    
    for (start in starts) {
        diffs=bisulf$pos-(start+2)<2
        motif=bisulf[diffs,]

        megadiffs=abs(mega$pos-start)<6
        megamotif=mega[megadiffs,]

        motifinfo=tibble(chr=chromname, pos=start, bisulf=mean(motif$bisulf), mega=max(megamotif$methfrac))
        methcompare=bind_rows(methcompare, motifinfo)
    }
}

    

boxpdf=file.path(dbxdir, 'bisulf_boxes.pdf')
pdf(boxpdf, h=9, w=13)
box=ggplot(methcompare, aes(x=chr, y=bisulf, colour=chr, fill=chr, alpha=.2)) +
    geom_boxplot() +
    ggtitle('GCCGGC') +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(box)
dev.off()

boxpdf=file.path(dbxdir, 'mega_boxes.pdf')
pdf(boxpdf, h=9, w=13)
box=ggplot(methcompare, aes(x=chr, y=mega, colour=chr, fill=chr, alpha=.2)) +
    geom_boxplot() +
    ggtitle('GCCGGC') +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(box)
dev.off()

corpdf=file.path(dbxdir, 'bisulf_cor.pdf')
pdf(corpdf, h=9, w=13)
box=ggplot(methcompare, aes(x=bisulf, y=mega, colour=chr, alpha=.2)) +
    geom_point() +
    ggtitle('GCCGGC') +
    theme_bw()
print(box)
dev.off()

boxpdf=file.path(dbxdir, 'mega_boxes_fixed.pdf')
pdf(boxpdf, h=9, w=13)
box=ggplot(methcomparefixed, aes(x=chr, y=mega, colour=chr, fill=chr, alpha=.2)) +
    geom_boxplot() +
    ggtitle('GCCGGC') +
    theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(box)
dev.off()

corpdf=file.path(dbxdir, 'bisulf_cor_fixed.pdf')
pdf(corpdf, h=9, w=13)
box=ggplot(methcomparefixed, aes(x=bisulf, y=mega, colour=chr, alpha=.2)) +
    geom_point() +
    ggtitle('GCCGGC') +
    theme_bw()
print(box)
dev.off()





##sanity check
allmethecoli=allmeth %>%
    filter(chr=='Escherichia_coli_chromosome')
ecolicg=meth %>%
    filter(chrom=='Escherichia_coli_chromosome') %>%
    filter(motif=='GCCGGC')

