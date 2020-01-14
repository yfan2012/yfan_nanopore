library(tidyverse)

##plot bulk qual scores for r9, r10

##plot qual score difference dists at places of overlap, no overlap

datadir='/uru/Data/Nanopore/projects/nivar/'
r9qcsv=paste0(datadir, 'error_quals/r9_quals.csv')
r10qcsv=paste0(datadir, 'error_quals/r10_quals.csv')

r9quals=read_csv(r9qcsv) %>%
    mutate(pore='r9')
r10quals=read_csv(r10qcsv) %>%
    mutate(pore='r10')

