library(tidyverse)
source('prolificans_asm_functions.R')

datadir='/pym/Data/Nanopore/projects/prolificans'
dbxdir='~/Dropbox/timplab_data/prolificans'


strains=c('st31', 'st90853', 'st5317')

##add busco info to asm summary
asmsumfile=file.path(dbxdir, 'qc', 'asmstats_final.csv')
asmsumcols=c('asmname','numtigs', 'n50', 'maxlen', 'minlen', 'total' )
asmsum=read_csv(asmsumfile, col_names=asmsumcols, col_types=cols()) %>%
    rowwise() %>%
    mutate(asm=sub('\\.final$', '', asmname))


buscoinfo=tibble(asm=as.character(),
                 single=as.integer(),
                 dup=as.integer(),
                 frag=as.integer(),
                 missing=as.integer())
busco_cols=c('buscoid', 'status', 'contig', 'start', 'end', 'score', 'length', 'url', 'description')
for (i in strains){
    buscodir=file.path(datadir, i, 'busco')
    buscoruns=basename(list.dirs(buscodir, recursive=FALSE))
    for (run in buscoruns) {
        buscotablefile=file.path(datadir, i,'busco',  run, run, 'run_sordariomycetes_odb10', 'full_table.tsv')
        sum=read_tsv(buscotablefile, comment='#', col_names=busco_cols, col_types=cols()) %>%
            summarise(single=sum(status=='Complete'),
                   missing=sum(status=='Missing'),
                   dup=sum(status=='Duplicated'),
                   frag=sum(status=='Fragmented')) %>%
            mutate(asm=run)
        buscoinfo=bind_rows(buscoinfo, sum)
    }
}
        
fullasmsum=asmsum %>%
    full_join(buscoinfo, by=c('asm'))
    
fullsumfile=file.path(dbxdir, 'qc', 'fullasminfo.csv')
write_csv(fullasmsum, fullsumfile)
