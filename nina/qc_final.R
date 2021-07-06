library(tidyverse)
library(Biostrings)
library(RColorBrewer)
source('prolificans_asm_functions.R')
datadir='/pym/Data/Nanopore/projects/prolificans'
dbxdir='~/gdrive/prolificans'
strains=c('st31', 'st90853', 'st5317')    



##checking average coverage, and excluding tigs that have very low avg cov to get final genomes
##write out final genomes
alltiginfo=tibble(tigname=as.character(),
                  length=as.integer(),
                  npavg=as.numeric(),
                  illavg=as.numeric(),
                  npnone=as.integer(),
                  illnone=as.integer(),
                  asm=as.character())

for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes_mitotrim')
    newgenomedir=file.path(datadir, i, 'genomes_final')
    prefixes=sub('\\.fasta$', '', list.files(genomedir, '.fasta$'))

    system(paste0('mkdir -p ', newgenomedir))
    
    for (prefix in prefixes) {
        npcovfile=file.path(datadir, i, 'cov', paste0(prefix, '.cov'))
        illcovfile=file.path(datadir, i, 'cov', paste0(prefix, '.illumina.cov'))
        cov=get_asm_cov(npcovfile, illcovfile)
        tiginfo=get_tigs_info(cov)

        npmed=median(tiginfo$npavg)
        illmed=median(tiginfo$illavg)
        
        keeptigs=tiginfo %>%
            filter(npavg>npmed*.9 & illavg>illmed*.9)

        asmfile=file.path(genomedir, paste0(prefix, '.fasta'))
        asm=readDNAStringSet(asmfile)

        newfix=paste0(str_split(prefix, '\\.')[[1]][1], '.', str_split(prefix, '\\.')[[1]][2])     
        newasm=asm[names(asm)[names(asm) %in% keeptigs$tigname]]
        newasmfile=file.path(datadir, i, 'genomes_final', paste0(newfix, '.final.fasta'))
        writeXStringSet(newasm, newasmfile, format='fasta')

        alltiginfo=bind_rows(alltiginfo, keeptigs)
    }
}




##check for telomeres
fwdtelo='AGGGTTAGGGTT'
revtelo='AACCCTAACCCT'
##single telo repeat looks like it's too short to use with just dumb matching.
##looking for repeat pairs, then multiply count later.

alltelos=tibble(chr=as.character(),
                length=as.numeric(),
                fwd=as.numeric(),
                rev=as.numeric(),
                asm=as.character(),
                fwdlow=as.integer(),
                fwdhigh=as.integer(),
                revlow=as.integer(),
                revhigh=as.integer())


for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes_final')
    prefixes=sub('\\.fasta$', '', list.files(genomedir, '.fasta$'))
    for (prefix in prefixes) {
        asmfile=file.path(genomedir, paste0(prefix, '.fasta'))
        telos=telocheck(asmfile, fwdtelo, revtelo) %>%
            mutate(fwd=fwd*2) %>%
            mutate(fwdlow=fwdlow*2) %>%
            mutate(fwdhigh=fwdhigh*2) %>%
            mutate(rev=rev*2) %>%
            mutate(revlow=revlow*2) %>%
            mutate(revhigh=revhigh*2) %>%
            mutate(asm=prefix)
        alltelos=bind_rows(alltelos, telos)
    }
}

teloinfofile=file.path(dbxdir, 'telo_info_final.csv')
write_csv(alltelos, teloinfofile)

fulltiginfo=inner_join(alltelos, alltiginfo, by=c('asm', 'length','chr'='tigname'))
tiginfofile=file.path(dbxdir, 'contig_info_final.csv')
write_csv(fulltiginfo, tiginfofile)

allteloplots=tibble(telopos=as.integer(),
                    percpos=as.numeric(),
                    telo=as.character(),
                    chr=as.character(),
                    name=as.character())

for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes_final')
    prefixes=sub('\\.fasta$', '', list.files(genomedir, '.fasta$'))
    for (prefix in prefixes) {
        asmfile=file.path(genomedir, paste0(prefix, '.fasta'))
        teloplotdf=as_tibble(teloplot(asmfile, fwdtelo, revtelo)) %>%
            mutate(name=prefix)

        allteloplots=bind_rows(teloplotdf, allteloplots)
    }
}
allteloplots=allteloplots %>%
    rowwise() %>%
    mutate(strain=str_split(name, '\\.')[[1]][1])

teloplotfile=file.path(dbxdir, 'telos_final.pdf')
pdf(teloplotfile, w=11, h=8)
for (i in strains) {
    straintelo=allteloplots %>%
        filter(strain==i)
    
    print(ggplot(straintelo, aes(x=percpos, colour=name, fill=name, alpha=.2)) +
        geom_histogram(position='identity') +
        scale_fill_brewer(palette = "Set2") +
        scale_colour_brewer(palette = "Set2") +
        facet_wrap(. ~ name, ncol=2) +
        ggtitle('Telo positions') +
        xlab('position (as percent of seq length)') +
        theme_bw())
}
dev.off()



