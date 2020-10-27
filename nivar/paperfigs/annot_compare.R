library(tidyverse)
##many bits referenced from gmoney's moth code. ty, gmoney.

datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/annotation'

gff=tibble(
    liftcer=file.path(datadir, 'liftoff', 'nivar_cer_lifted.gff'),
    ##liftalb=file.path(datadir, 'liftoff', 'nivar_alb_lifted.gff'),
    liftgla=file.path(datadir, 'liftoff', 'nivar_gla_lifted.gff'),
    braker=file.path(datadir, 'braker', 'braker.gff3'),
    drna=file.path(datadir, 'stringtie', 'denovo_drna.gff'),
    rnaseq=file.path(datadir, 'stringtie', 'denovo_rnaseq.gff'))


cols=c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
liftcergff=read_tsv(gff$liftcer, col_names=cols)
##liftalbgff=read_tsv(gff$liftalb, col_names=cols)
liftglagff=read_tsv(gff$liftgla, col_names=cols)
brakergff=read_tsv(gff$braker, col_names=cols)
drnagff=read_tsv(gff$drna, col_names=cols, skip=2)


##find annots included in liftoff albicans that weren't found in liftoff cerevisiae
'''
liftalbfile=file.path(datadir, 'liftoff', 'liftcer_liftalb.nivar_alb_lifted.gff.tmap')
liftalb=read_tsv(liftalbfile) %>%
    filter(class_code=='u')
fromalb=liftalbgff %>%
    rowwise() %>%
    mutate(geneid=substring(strsplit(attribute, ';', fixed=TRUE)[[1]][1], 4)) %>%
    filter(geneid %in% liftalb$qry_gene_id | geneid %in% liftalb$qry_id) %>%
    select(-geneid)
lifted_annots=rbind(liftcergff, fromalb)

lifted_tsv=file.path(datadir, 'combined', 'lifted_all.gff')
write_tsv(lifted_annots, lifted_tsv, col_names=FALSE)
'''

liftcerfile=file.path(datadir, 'liftoff', 'liftgla_liftcer.nivar_cer_lifted.gff.tmap')
liftcer=read_tsv(liftcerfile) %>%
    filter(class_code=='u')
fromcer=liftcergff %>%
    rowwise() %>%
    mutate(geneid=substring(strsplit(attribute, ';', fixed=TRUE)[[1]][1], 4)) %>%
    filter(geneid %in% liftcer$qry_gene_id | geneid %in% liftgla$qry_id) %>%
    select(-geneid)
lifted_annots=rbind(liftcergff, fromgla)

##lifted_tsv=file.path(datadir, 'combined', 'lifted_all.gff')
lifted_tsv=file.path(datadir, 'combined', 'lifted_all_gla.gff')
write_tsv(lifted_annots, lifted_tsv, col_names=FALSE)

##find common ones between stringtie and braker
drnabrakerfile=file.path(datadir, 'braker', 'drna_braker.braker.gff3.tmap')
drnabraker=read_tsv(drnabrakerfile) %>% #keep these from drna annot
    filter(class_code=='=' | class_code=='c')
drnabraker_rev=read_tsv(drnabrakerfile) %>% #keep these from braker annot
    filter(class_code=='k')

fromdrna=drnagff %>%
    rowwise() %>%
    mutate(geneid=strsplit(attribute, '"', fixed=TRUE)[[1]][2]) %>%
    filter(geneid %in% drnabraker$ref_gene_id | geneid %in% drnabraker$ref_id) %>%
    select(-geneid)
frombraker=brakergff %>%
    rowwise() %>%
    mutate(geneidraw=substring(strsplit(attribute, ';', fixed=TRUE)[[1]][1], 4)) %>%
    mutate(geneid=paste0(strsplit(geneidraw, '.', fixed=TRUE)[[1]][1],'.', strsplit(geneidraw, '.', fixed=TRUE)[[1]][2])) %>%
    filter(geneid %in% drnabraker_rev$qry_gene_id | geneid %in% drnabraker_rev$qry_id) %>%
    #mutate(attribute=paste0('"', attribute, '"')) %>%
    select(-geneid, -geneidraw)
data_annots=rbind(frombraker, fromdrna)

data_tsv=file.path(datadir, 'combined', 'data_all.gff')
write.table(data_annots, data_tsv, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)


