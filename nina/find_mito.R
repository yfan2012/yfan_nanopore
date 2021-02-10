library(tidyverse)
library(Biostrings)
library(RColorBrewer)
library(BSgenome)
library(GenomicRanges)
source('prolificans_asm_functions.R')

datadir='/pym/Data/Nanopore/projects/prolificans'
dbxdir='~/Dropbox/timplab_data/prolificans'

strains=c('st31', 'st90853', 'st5317')

teloinfofile=file.path(dbxdir, 'telo_info.csv')
alltelos=read_csv(teloinfofile, col_types=cols())

##trim mitos
for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes_covfilt')
    prefixes=sub('\\.covfilt.fasta$', '', list.files(genomedir, '.covfilt.fasta$'))

    newgenomedir=file.path(datadir, i, 'genomes_mitotrim')
    system(paste0('mkdir -p ', newgenomedir))

    
    for (asm in prefixes) {
        coordsfile=file.path(datadir, i, 'mummer_mito', paste0(asm, '.mcoords'))
        asmtiginfo=alltelos %>%
            filter(asm==asm)
        asmfile=file.path(genomedir, paste0(asm, '.covfilt.fasta'))
        outfile=file.path(newgenomedir, paste0(asm, '.mitotrim.fasta'))
        print(asm)
        suggest_mito_breaks(coordsfile, asmtiginfo, asmfile, outfile)
    }
}
