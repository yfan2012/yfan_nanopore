library(tidyverse)
library(Biostrings)

datadir='/pym/Data/Nanopore/projects/prolificans/st31'
finalfile=file.path(datadir, 'final/st31.final.fasta')
genome=readDNAStringSet(finalfile, format='fasta')

####put back tig 20 and tig 13 because repeat regions there explain low coverage
unbrokenfile=file.path(datadir, 'genomes/st31.ragtag_fc.fasta')
unbroken=readDNAStringSet(unbrokenfile, format='fasta')

keep=genome[-c(2,3,8,9)]
get=unbroken[c(2,7)]
newgenome=c(keep, get)
newgenomefile=file.path(datadir, 'final/st31.final2.fasta')
writeXStringSet(newgenome, newgenomefile)
