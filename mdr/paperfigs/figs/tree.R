library(tidyverse)
library(multidplyr)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(fossil)
library(mclust)
library(ggtree)
library(treeio)
library(ape)
source('~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clinical_functions.R')

projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')

dbxdir='~/gdrive/mdr/paperfigs/figs'


##chrom to bin info
chrombinsfile='/mithril/Data/Nanopore/projects/methbin/paperfigs/contig_level/tigs2bins.tsv'
chrombins=read_tsv(chrombinsfile)
##chrom to taxnonomy info
taxfile=file.path('~/gdrive/mdr/paperfigs/contig_level/tigbins_species.csv')
taxinfo=read_csv(taxfile)

####get tree from methylation distance
exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG', 'GCCGGC', 'TCCGGA')
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))
fullfreqs=methfreq %>%
    rowwise() %>%
    spread(key=motif, value=freq)
fullnacount=colSums(is.na(fullfreqs))/dim(fullfreqs)[1]
fullcounts=tibble(motifs=names(fullnacount), frac=1-fullnacount)

fullcountcsv=file.path(dbxdir, 'clin_motifcounts.csv')
write_csv(fullcounts, fullcountcsv)


freqs=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec) %>%
    spread(key=motif, value=freq)
nacount=colSums(is.na(freqs))/dim(freqs)[1]

nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])

methchroms=methfreq %>%
    filter(chrom %in% keepchroms)
chrominfo=methchroms %>%
    spread(key=motif, value=freq) %>%
    filter(chrom %in% chrombins$rname)
matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom
dend=matchrominfo %>%
    scale %>% 
    dist %>%
    hclust %>%
    as.dendrogram
dendfile=file.path(dbxdir, 'plaindend.tree')
write.tree(dend %>% as.phylo, dendfile)


labelinfo=tibble(label=labels(dend)) %>%
    filter(label %in% chrombins$rname) %>%
    rowwise() %>%
    mutate(bins=chrombins$bin[chrombins$rname==label])


mycolors=read_csv('colors.csv')

for (i in labelinfo$label) {
    if (!i %in% taxinfo$tig) {
        bin=labelinfo$bins[labelinfo$label==i]
        binleaf=taxinfo$binleaf[taxinfo$bin==bin][1]
        unkinfo=tibble(tig=i, bin=bin, tigleaf='Not classified', binleaf=binleaf, lcalevel='unk', lca='unk')
        taxinfo=bind_rows(taxinfo, unkinfo)
    }
}


labelinfo=labelinfo %>%
    mutate(color=mycolors$color[mycolors$bin==bins]) %>%
    mutate(tax=taxinfo$tigleaf[taxinfo$tig==label])

plaindend=dend %>%
    set('leaves_pch', 15) %>%
    set('leaves_cex', 2) %>%
    set('leaves_col', labelinfo$color) %>%
    set('labels', rep('', length(labelinfo$color)))
    ##set('labels', labelinfo$tax) %>%
    ##set('by_labels_branches_col', labelinfo$color)



treepdf=file.path(dbxdir, 'tree.pdf')
pdf(treepdf, h=35, w=8)
par(mar = c(3, 4, 2, 15) + 0.1)
plot(plaindend, horiz=TRUE)
dev.off()


##show plasmid labels
plasnearbinscsv=file.path(dbxdir, 'plasmid_nearest.csv')
plasnear=read_csv(plasnearbinscsv)
labelinfo=labelinfo %>%
    mutate(type=case_when(label %in% unique(plasnear$chroms) ~ label,
                          TRUE ~ ''))

dendlabeled=dend %>%
    set('leaves_pch', 15) %>%
    set('leaves_cex', 2) %>%
    set('leaves_col', labelinfo$color) %>%
    set('labels', labelinfo$type)

labtreepdf=file.path(dbxdir, 'tree_labeled.pdf')
pdf(labtreepdf, h=35, w=8)
par(mar = c(3, 4, 2, 15) + 0.1)
plot(dendlabeled, horiz=TRUE)
dev.off()
