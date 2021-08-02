library(tidyverse)
source('prolificans_asm_functions.R')
library(ggforce)

datadir='/pym/Data/Nanopore/projects/prolificans'
dbxdir='~/gdrive/prolificans'


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




###make busco plots
busco_cols=c('buscoid', 'status', 'contig', 'start', 'end', 'score', 'length', 'url', 'description')
allbuscos=tibble(
    status=as.character(),
    sum=as.numeric(),
    asm=as.character())
for (i in strains) {
    buscofile=file.path(datadir, i, 'final/busco/run_sordariomycetes_odb10/full_table.tsv')
    busco=read_tsv(buscofile, comment='#', col_names=busco_cols) %>%
        group_by(status) %>%
        summarise(sum=length(unique(buscoid))) %>%
        mutate(asm=i)
    allbuscos=bind_rows(allbuscos, busco)
}

buscoplot=file.path(dbxdir, 'busco_genome.pdf')
pdf(buscoplot, height=8, width=19)
barplot=ggplot(allbuscos, aes(x=sum, y=asm, fill=status, alpha=.8)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    ggtitle('BUSCO') +
    ylab('Genome') +
    xlab('Number of BUSCOs') +
    theme_bw()
print(barplot)
dev.off()

buscozoom=file.path(dbxdir, 'busco_genome_zoom.pdf')
pdf(buscozoom, height=8, width=15)
print(barplot+facet_zoom(xlim=c(0,300)))
dev.off()

tscriptbuscos=tibble(
    status=as.character(),
    sum=as.numeric(),
    asm=as.character())
transcriptomes=c('st5137_braker', 'st5137_trinity', 'st5137_transcriptome')
for (i in transcriptomes) {
    buscofile=file.path(datadir, 'rna/busco', i, 'run_sordariomycetes_odb10/full_table.tsv')
    busco=read_tsv(buscofile, comment='#', col_names=busco_cols) %>%
        group_by(status) %>%
        summarise(sum=length(unique(buscoid))) %>%
        mutate(asm=i)
    tscriptbuscos=bind_rows(tscriptbuscos, busco)
}

tscriptfile=file.path(dbxdir, 'busco_transcriptome.pdf')
pdf(tscriptfile, height=8, width=15)
tscriptplot=ggplot(tscriptbuscos, aes(x=sum, y=asm, fill=status, alpha=.8)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    ggtitle('BUSCO') +
    ylab('Genome') +
    xlab('Number of BUSCOs') +
    theme_bw()
print(tscriptplot)
dev.off()




##missing buscos
missing=NULL
for (i in strains) {
    buscofile=file.path(datadir, i, 'final/busco/run_sordariomycetes_odb10/full_table.tsv')
    busco=read_tsv(buscofile, comment='#', col_names=busco_cols) %>%
        mutate(asm=i) %>%
        filter(status=='Missing')
    missing=bind_rows(missing, busco)
}

library(jsonlite)
busco_description <- function(busco_id) {
    ##from https://thackl.github.io/BUSCO-gene-descriptions
    ##get busco info
    odb_info=read_json(paste0("https://www.orthodb.org/group?id=", busco_id), simplifyVector = TRUE)
    info=odb_info$data$name
    return(info)
}
        
allmissing=missing %>%
    group_by(buscoid) %>%
    filter(n()==3) %>%
    summarise(id=buscoid[1]) %>%
    select(-id) %>%
    rowwise() %>%
    mutate(description=busco_description(buscoid))

missingfile=file.path(dbxdir, 'missing_buscos_from_all.txt')
write_tsv(allmissing, missingfile)

tscriptmissing=NULL
for (i in transcriptomes) {
    buscofile=file.path(datadir, 'rna/busco', i, 'run_sordariomycetes_odb10/full_table.tsv')
    busco=read_tsv(buscofile, comment='#', col_names=busco_cols) %>%
        mutate(asm=i) %>%
        filter(status=='Missing')
    tscriptmissing=bind_rows(tscriptmissing, busco)
}
tscriptmissing=tscriptmissing %>%
    group_by(buscoid) %>%
    summarise(id=buscoid[1]) %>%
    select(-id) %>%
    rowwise() %>%
    mutate(description=busco_description(buscoid))
missingtscriptfile=file.path(dbxdir, 'missing_buscos_from_all_transcriptome.tsv')
write_tsv(tscriptmissing, missingtscriptfile)

i="st5137_trinity"
buscofile=file.path(datadir, 'rna/busco', i, 'run_sordariomycetes_odb10/full_table.tsv')
busco=read_tsv(buscofile, comment='#', col_names=busco_cols) %>%
    mutate(asm=i) %>%
    filter(status=='Missing')
trinitymissing=busco %>%
    select(buscoid) %>%
    rowwise() %>%
    mutate(description=busco_description(buscoid))
missingtrinityfile=file.path(dbxdir, 'missing_buscos_trinity.tsv')
write_tsv(trinitymissing, missingtrinityfile)

missingall=trinitymissing %>%
    rowwise() %>%
    filter(buscoid %in% allmissing$buscoid) %>%
    filter(buscoid %in% tscriptmissing$buscoid)
totalmissing=file.path(dbxdir, 'missing_buscos_genome_and_transcriptome.tsv')
write_tsv(missingall, totalmissing)




allbuscos=tibble(
    status=as.character(),
    sum=as.numeric(),
    asm=as.character())
for (i in strains) {
    buscofile=file.path(datadir, i, 'busco', i, i, '/run_sordariomycetes_odb10/full_table.tsv')
    busco=read_tsv(buscofile, comment='#', col_names=busco_cols) %>%
        group_by(status) %>%
        summarise(sum=length(unique(buscoid))) %>%
        mutate(asm=i)
    allbuscos=bind_rows(allbuscos, busco)
}

buscoplot=file.path(dbxdir, 'busco_genome_final.pdf')
pdf(buscoplot, height=8, width=19)
barplot=ggplot(allbuscos, aes(x=sum, y=asm, fill=status, alpha=.8)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    ggtitle('BUSCO') +
    ylab('Genome') +
    xlab('Number of BUSCOs') +
    theme_bw()
print(barplot)
dev.off()

buscofinalcsv=file.path(dbxdir, 'busco_genome_final.csv')
write_csv(allbuscos, buscofinalcsv)
