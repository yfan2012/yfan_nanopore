library(tidyverse)
library(Biostrings)
library(BSgenome)

datadir='/pym/Data/Nanopore/projects/prolificans/st5317'
finalfile=file.path(datadir, 'final/st5317.final.fasta')

##read current genome
genome=readDNAStringSet(finalfile, format='fasta')


####correct inserted mito in tig14
##read mummer of current genome
mitomumfile=file.path(datadir, 'mummer_mito/st5317_correction.mcoords')
mumheads=c('rstart', 'rend', 'qstart', 'quend', 'ralen', 'qalen', 'ident', 'rlen', 'alen', 'sim', 'cov', 'rname', 'qname')
mum=read_tsv(mitomumfile, col_names=mumheads) %>%
    filter(rname=='tig14_RagTag')

##take out last 6kb of tig14. confirmed coverage jump in IGV
newtig=DNAStringSet(substr(as.character(genome['tig14_RagTag']), 1, mum$rstart[1]-1))
names(newtig)=c('tig14_RagTag_corrected')


####get full length mito
##transplant canu mito into final
canufile=file.path(datadir, 'genomes_mitotrim/st5317.canu.mitotrim.fasta')
canu=readDNAStringSet(canufile, format='fasta')
newmito=canu['st5317.canu_mito']


####get missing section
untrimmedfile=file.path(datadir, 'genomes_mitotrim/st5317.ragtag_cf.mitotrim.fasta')
untrimmed=readDNAStringSet(untrimmedfile, format='fasta')
missing=untrimmed["tig00000038_RagTag_break1"]
names(missing)=c("tig38_RagTag_break1")


####cat together
newgenome=c(genome[-c(7,17)], newmito, missing)
newgenomefile=file.path(datadir, 'final/st5317.final2.fasta')
writeXStringSet(newgenome, newgenomefile)

