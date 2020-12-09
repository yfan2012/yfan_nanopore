library(tidyverse)
library(tidyr)
library(gridExtra)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(R.utils)
library(foreach)
library(doParallel)

dbxdir='~/Dropbox/yfan/nivar/'
datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/'

reffile='/uru/Data/Nanopore/projects/nivar/reference/candida_nivariensis.fa'
glafile='/uru/Data/Nanopore/projects/nivar/reference/medusa_fungi/candida_glabrata.fa'
xufile='/uru/Data/Nanopore/projects/nivar/paperfigs/patho/glabrata_xu.fa'

regionsfile='/uru/Data/Nanopore/projects/nivar/paperfigs/patho/gpicwp_regions.csv'


##read in regions file from xu (asterisks removed)
regiontable=read_csv(regionsfile) %>%
    filter(!str_detect(Status, 'Removed')) %>%
    select(-Source, -Note, -CommonName)
regions=GRanges(seqnames=regiontable$Chrom,
                ranges=IRanges(start=regiontable$Start, end=regiontable$End),
                strand=regiontable$Direction,
                gene=regiontable$Systematic_Name)

##rename so the chrs match the table
xugla=readDNAStringSet(xufile)
chrnames=tibble(chr='Chr', letter=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M')) %>%
    mutate(chrname=paste0(chr, letter))
names(xugla)=chrnames$chrname

##get gene seqs
geneseqs=getSeq(xugla, regions)
names(geneseqs)=regiontable$Systematic_Name

##write to fasta
geneseqsfa=file.path(datadir, 'patho', 'gpicwp.fa')
writeXStringSet(geneseqs, geneseqsfa, format='fasta')

