library(tidyverse)
library(Biostrings)
library(UpSetR)

datadir='/uru/Data/Nanopore/projects/mdr'
dbxdir='~/Dropbox/timplab_data/mdr'
rebasefiltcsv=file.path(datadir, 'refs', 'rebase_filtered.csv')
zymoreffa='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'

rebase=read_csv(rebasefiltcsv)
zymoref=readDNAStringSet(zymoreffa)
species=c('Pseudomonas aeruginosa', 'Escherichia coli' )


##hand curating the strains, becuase they're listed somewhat inconsistently
bs=rebase %>%
    rowwise() %>%
    filter(grepl('Bacillus subtilis', Organism, fixed=TRUE)) %>%
    filter(Strain=='ATCC 6633')
##none found for this strain. can't find other names either
ef=rebase %>%
    rowwise() %>%
    filter(grepl('faecalis', Organism, fixed=TRUE)) #%>%
    #filter(Strain=='ATCC 7080')
##https://ecoliwiki.org/colipedia/index.php/Escherichia_coli_(Migula_1895)_Castellani_and_Chalmers_1919
ec=rebase %>%
    rowwise() %>%
    filter(grepl('Escherichia coli', Organism, fixed=TRUE)) %>%
    filter(Strain %in% c('DSM 30083 = JCM 1649 = ATCC 11775','NCTC9001'))
lm=rebase %>%
    rowwise() %>%
    filter(grepl('Listeria monocytogenes', Organism, fixed=TRUE)) %>%
    filter(Strain=='ATCC 19117')
pa=rebase %>%
    rowwise() %>%
    filter(grepl('Pseudomonas aeruginosa', Organism, fixed=TRUE)) %>%
    filter(Strain=='ATCC 15442')
se=rebase %>%
    rowwise() %>%
    filter(grepl('Salmonella enterica', Organism, fixed=TRUE)) %>%
    filter(grepl('13311', Strain, fixed=TRUE))
##looking up this one with the dsm number
sa=rebase %>%
    rowwise() %>%
    filter(grepl('Staphylococcus aureus', Organism, fixed=TRUE)) %>%
    filter(grepl('20231', Strain, fixed=TRUE))


zymo=bind_rows(bs,ec,lm,pa,se,sa) %>%
    filter(Type=='II', nchar(Specificity)>3)
motifs=unique(zymo$Specificity)


library(BMS)
barcode_num <- function(org, motifs, zymo) {
    orgzymo=zymo %>%
        filter(Organism==org)
    barcode=motifs %in% orgzymo$Specificity
    bc=bin2hex(as.integer(barcode))
}

bcinfo=tibble(org=unique(zymo$Organism)) %>%
    rowwise() %>%
    mutate(bc=barcode_num(org, motifs, zymo))



rebaseii=rebase %>%
    filter(Type=='II', nchar(Specificity)>3)
allmotifs=unique(rebaseii$Specificity)
infoii=tibble(org=unique(rebaseii$Organism)) %>%
    rowwise() %>%
    mutate(bc=barcode_num(org, allmotifs, rebaseii)) %>%
    group_by(bc) %>%
    summarise(num=n()) %>%
    arrange(-num)


    









