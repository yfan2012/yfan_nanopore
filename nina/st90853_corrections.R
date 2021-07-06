library(tidyverse)
library(Biostrings)

datadir='/pym/Data/Nanopore/projects/prolificans/st90853'
finalfile=file.path(datadir, 'final/st90853.final.fasta')
genome=readDNAStringSet(finalfile, format='fasta')

####get rid of tigs 16,17,18 because they're redundant with each other and scaffold19
newgenome=genome[-c(12,13,14)]
newgenomefile=file.path(datadir, 'final/st90853.final2.fasta')
writeXStringSet(newgenome, newgenomefile)
