library(tidyverse)

datadir='~/data/mdr/zymo/barcode'
dbxdir='~/Dropbox/timplab_data/mdr/zymo'

barcodefile=file.path(datadir, '20190809_zymo_control_barcodes_filtered.txt')
bc_cols=c('readname', 'chrname', 'start', 'end', 'mapq', 'GATC', 'CCWGG', 'ATGCAT','GTCGAC','CTCCAG','CTKVAG')

bcinfo=read_tsv(barcodefile, col_names=bc_cols) %>%
    select(-start, -end, -mapq)


