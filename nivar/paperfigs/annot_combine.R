library(tidyverse)
##many bits referenced from gmoney's moth code. ty, gmoney.
rawdir='/uru/Data/Nanopore/projects/nivar/paperfigs'
datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/annotation'

cols=c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
lifted_tsv=file.path(datadir, 'combined', 'lifted_all.gff')
data_tsv=file.path(datadir, 'combined', 'data_all.gff')


##gffcompare between lifted_annots and data_annots
deltfile=file.path(datadir, 'combined', 'lifted_vs_data.data_all.gff.tmap')
delt=read_tsv(deltfile) %>%
    filter(class_code=='u')

lifted=read_tsv(lifted_tsv, col_names=cols)
brakerdata=read_tsv(data_tsv, col_names=cols) %>%
    filter(source=='AUGUSTUS') %>%
    rowwise() %>%
    mutate(geneidraw=substring(strsplit(attribute, ';', fixed=TRUE)[[1]][1], 4)) %>%
    mutate(geneid=paste0(strsplit(geneidraw, '.', fixed=TRUE)[[1]][1],'.', strsplit(geneidraw, '.', fixed=TRUE)[[1]][2])) %>%
    filter(geneid %in% delt$qry_gene_id | geneid %in% delt$qry_id) %>%
    select(-geneid, -geneidraw)
drnadata=read_tsv(data_tsv, col_names=cols) %>%
    filter(source=='StringTie') %>%
    rowwise() %>%
    mutate(geneid=strsplit(attribute, '"', fixed=TRUE)[[1]][2]) %>%
    filter(geneid %in% delt$qry_gene_id | geneid %in% delt$qry_id) %>%
    select(-geneid)

final_annot=rbind(lifted, brakerdata, drnadata)
final_annot_tsv=file.path(rawdir,'annotation_final', 'nivar.final.gff')
write.table(final_annot, final_annot_tsv, row.names=FALSE, col.names=FALSE, sep='\t')
