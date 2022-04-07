library(Biostrings)
library(tidyverse)

reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
ref=readDNAStringSet(reffile)[c(1:6, 58:62)]
labelchrs=c('bsubtilis', 'efaecalis',rep('ecoli',2), 'lmonocytogenes', 'paeruginosa', 'senterica', rep('saureus', 4))
chrlabs=tibble(chrname=names(ref), label=labelchrs)

##grab bisulf data
dbxdir='~/gdrive/mdr/paperfigs/figs'
projdir='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite'
datadir=file.path(projdir, 'bismark')
labels=c('bsubtilis', 'efaecalis','ecoli', 'lmonocytogenes', 'nasa_Ecoli_K12', 'paeruginosa', 'saureus', 'senterica')
cxcols=c('chr', 'pos', 'strand', 'meth', 'unmeth', 'context', 'seq')

##motif info
motiffile=file.path(projdir, 'zymo_cmeth_expanded.csv')
motifinfo=read_csv(motiffile)



##dealing with chr and plasmid separately
chromcounts=NULL
for (i in names(ref)) {
    label=chrlabs$label[chrlabs$chrname==i]
    
    cxfile=file.path(datadir, label, paste0(label, '_1_bismark_bt2_pe.CX_report.txt'))
    cx=read_tsv(cxfile, col_names=cxcols) %>%
        filter(!grepl('tig', chr, fixed=TRUE)) %>%
        filter((meth+unmeth)>15) %>%
        mutate(methfrac=meth/(meth+unmeth))

    motifs=motifinfo$motif[motifinfo$species==label]

    shorti=strsplit(i, ' ')[[1]][1]
    cxchr=cx %>%
        filter(chr==shorti)
    
    methylated=cxchr %>%
        filter(methfrac>.9)

    methknown=NULL
    for (motif in motifs) {
        chromseq=ref[i]
        motifmatches=vmatchPattern(motif, chromseq)[[1]]

        methmatches=methylated %>%
            rowwise() %>%
            mutate(inmotif=case_when(length(findOverlaps(IRanges(start=pos, end=pos), motifmatches))>0 ~ 'yes',
                                     TRUE ~ 'no')) %>%
            filter(inmotif=='yes') %>%
            mutate(motif=motif)
        methknown=bind_rows(methknown, methmatches)
    }

    motifcount=table(methknown$motif)
    
    nameinfo=tibble(chrom=i,
                    motifs=c(names(motifcount), 'none'),
                    counts=c(motifcount, dim(methylated)[1]-dim(methknown)[1]))

    chromcounts=bind_rows(chromcounts, nameinfo)
}                    

    
outfile=file.path(dbxdir, 'bisulfite_motif_counts.csv')
write_csv(chromcounts, outfile)
