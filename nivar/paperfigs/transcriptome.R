library(tidyverse)
library(tidyr)
library(gridExtra)
library(Biostrings)
library(R.utils)
library(foreach)
library(doParallel)

dbxdir='~/Dropbox/yfan/nivar/paperfigs/raw'
datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/'
finfile=paste0(datadir,'assembly_final/nivar.final.fasta')
finfilegff=paste0(datadir, 'annotation_final/nivar.final.gff')

reffile='/uru/Data/Nanopore/projects/nivar/reference/candida_nivariensis.fa'

glafile='/uru/Data/Nanopore/projects/nivar/reference/candida_glabrata.fa'
glafilegff='/uru/Data/Nanopore/projects/nivar/reference/candida_glabrata.gff'

albfile='/uru/Data/Nanopore/projects/nivar/reference/candida_albicans.fa'
albfilegff='/uru/Data/Nanopore/projects/nivar/reference/candida_albicans.gff'

cerfile='/uru/Data/Nanopore/projects/nivar/reference/saccharomyces_cerevisiae.fa'
cerfilegff='/uru/Data/Nanopore/projects/nivar/reference/saccharomyces_cerevisiae.gff'

cnames=c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')


all=read_tsv(nivargff, col_names=cnames, comment='#')
table(all$source)

get_ref_exons <- function(gfffile) {
    ##check for multiexonic genes
    gff=read_tsv(gfffile, col_names=cnames, comment='#') %>%
        filter(feature=='exon') %>%
        separate(attribute, into=c('prelocus', 'postlocus'), sep='locus_tag=', extra='drop') %>%
        separate(postlocus, into=c('locus'), sep=';', extra='drop')
    numexons=dim(gff)[1]
    numgenes=length(unique(gff$locus))
    return(c(numexons, numgenes))
}

get_nivar_exons  <- function(nivargff) {
    ##check for multiexonic genes
    ##different text parsing for AUGUSTUS, GeneMark.hmm, Liftoff, StringTie
    lift=read_tsv(nivargff, col_names=cnames, comment='#') %>%
        filter(source=='Liftoff') %>%
        filter(feature=='exon') %>%
        separate(attribute, into=c('prelocus', 'postlocus'), sep='locus_tag=', extra='drop') %>%
        separate(postlocus, into=c('gene'), sep=';', extra='drop')
    aug=read_tsv(nivargff, col_names=cnames, comment='#') %>%
        filter(source=='AUGUSTUS') %>%
        filter(feature=='exon') %>%
        separate(attribute, into=c('gene'), sep='\\.', extra='drop') %>%
        rowwise() %>%
        mutate(gene=substring(gene, 4))
    gm=read_tsv(nivargff, col_names=cnames, comment='#') %>%
        filter(source=='GeneMark.hmm') %>%
        filter(feature=='exon') %>%
        separate(attribute, into=c('gene'), sep='\\.', extra='drop') %>%
        rowwise() %>%
        mutate(gene=substring(gene, 4))
    st=read_tsv(nivargff, col_names=cnames, comment='#') %>%
        filter(source=='StringTie') %>%
        filter(feature=='exon') %>%
        separate(attribute, into=c('geneid'), sep=';', extra='drop') %>%
        separate(geneid, into=c('extra', 'gene'), sep='\\s', extra='drop') %>%
        rowwise() %>%
        mutate(gene=noquote(gene))
    exons=sum(dim(lift)[1], dim(aug)[1], dim(gm)[1], dim(st)[1])
    genes=sum(length(unique(lift$gene)), length(unique(aug$gene)), length(unique(gm$gene)), length(unique(st$gene)))
    return(c(exons, genes))
}

glacounts=get_ref_exons(glafilegff)
cercounts=get_ref_exons(cerfilegff)
albcounts=get_ref_exons(albfilegff)
fincounts=get_nivar_exons(finfilegff)


                                            
