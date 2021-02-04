library(tidyverse)
library(Biostrings)
library(RColorBrewer)
datadir='/pym/Data/Nanopore/projects/prolificans'
dbxdir='~/Dropbox/timplab_data/prolificans'

##give some contig stats
get_asm_cov <- function(npcovfile, illcovfile) {
    ##get coverage info for the asm
    filename=str_split(basename(npcovfile), '\\.')[[1]]
    asm=paste0(filename[1:length(filename)-1], collapse='.')
        
    covnames=c('tigname', 'start', 'end', 'pos', 'cov')
    npcov=read_tsv(npcovfile, col_names=covnames) %>%
        select(-start) %>%
        mutate(reads='np') %>%
        mutate(asm=asm)
    illcov=read_tsv(illcovfile, col_names=covnames) %>%
        select(-start) %>%
        mutate(reads='ill') %>%
        mutate(asm=asm)

    cov=rbind(npcov, illcov)
    return(cov)
}


get_tigs_info <- function(cov) {
    ##make table of info per contig
    ##tigname, length, mean np cov, mean illumina cov, # of no cov illumina, # of no cov nanopore
    ##number of likely breaks?? separate function
    tiginfo=cov %>%
        group_by(tigname) %>%
        summarise(length=mean(end),
                  npavg=mean(cov[reads=='np']),
                  illavg=mean(cov[reads=='ill']),
                  npnone=sum(cov[reads=='np']==0),
                  illnone=sum(cov[reads=='ill']==0),
                  asm=asm[1])
    return(tiginfo)
}


findbreaks <- function(tigcov) {
    ##given a chunk of zero cov data frame for , return tigbreaks
    ##get positions that are zero cov in both illumina and nanopore
    breaks=tigcov %>%
        group_by(pos) %>%
        summarise(tigname=tigname[1],
                  pos=pos[1],
                  support=n()) %>%
        filter(support==2)
    
    ##collapse into ranges
    diffs=diff(breaks$pos)

    start=c(breaks$pos[1])
    end=c()
    breaklocs=which(diffs!=1)

    for (i in breaklocs) {
        end=c(end, breaks[i,]$pos)
        start=c(start, breaks[i+1,]$pos)
    }
    end=c(end, breaks$pos[dim(breaks)[1]])

    breakranges=tibble(start=start, end=end)
    
    return(breakranges)
}

suggest_breaks <- function(cov) {
    ##look at spots of zero cov, and figure out breaks
    nocov=cov %>%
        filter(cov==0) %>%
        group_by(asm, tigname) %>%
        do(findbreaks(.))

    return(nocov)
}

telocheck <- function(asmfile, fwdtelo, revtelo) {
    ##snatched from nivar analysis
    asm=readDNAStringSet(asmfile)
    seqs=as.character(asm)
    
    teloinfo=tibble(chr=names(asm)) %>%
        mutate(length=width(asm[chr])) %>%
        mutate(fwd=str_count(as.character(asm[chr]), fwdtelo)) %>%
        mutate(rev=str_count(as.character(asm[chr]), revtelo))
}

library(doParallel)
teloplot <- function(asmfile, fwdtelo, revtelo) {
    ##also stolen from nivar analysis
    asm=readDNAStringSet(asmfile)

    teloinfo=foreach(i=1:length(asm), .combine=rbind) %dopar% {
        seq=as.character(asm[i])
        len=width(asm[i])
        
        fwdlocs=gregexpr(fwdtelo, seq)[[1]]
        fwdinfo=data.frame(telopos=fwdlocs, percpos=fwdlocs/len, telo=as.character('fwd'), chr=as.character(names(asm[i])))
        revlocs=gregexpr(revtelo, seq)[[1]]
        revinfo=data.frame(telopos=revlocs, percpos=revlocs/len, telo=as.character('rev'), chr=as.character(names(asm[i])))

        allinfo=rbind(fwdinfo, revinfo)
        return(allinfo)
    }
    cleanteloinfo=as_tibble(teloinfo) %>%
        mutate(name=asmfile) %>%
        mutate(telo=as.character(telo)) %>%
        mutate(chr=as.character(chr)) %>%
        filter(telopos!=-1)
    return(cleanteloinfo)
}   


##checking coverage
strains=c('st31', 'st90853', 'st5317')
alltiginfo=tibble(tigname=as.character(),
                  length=as.integer(),
                  npavg=as.numeric(),
                  illavg=as.numeric(),
                  npnone=as.integer(),
                  illnone=as.integer(),
                  asm=as.character())
allbreaks=tibble(asm=as.character(),
                 tigname=as.character(),
                 start=as.numeric(),
                 end=as.numeric())

for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes')
    prefixes=sub('\\.fasta$', '', list.files(genomedir, '.fasta$'))
    
    for (prefix in prefixes) {
        npcovfile=file.path(datadir, i, 'cov', paste0(prefix, '.cov'))
        illcovfile=file.path(datadir, i, 'cov', paste0(prefix, '.illumina.cov'))
        
        cov=get_asm_cov(npcovfile, illcovfile)
        tiginfo=get_tigs_info(cov)
        breakregions=suggest_breaks(cov)

        if (dim(breakregions)[1]<1) {
            breakregions=tibble(asm=tiginfo$asm[1], tigname='none', start=0, end=0)
        }
        
        alltiginfo=bind_rows(alltiginfo, tiginfo)
        allbreaks=bind_rows(allbreaks, breakregions)
    }
}

breakinfofile=file.path(dbxdir, 'zero_cov.csv')
write_csv(allbreaks,breakinfofile)



##check for telomeres
fwdtelo='AGGGTTAGGGTT'
revtelo='AACCCTAACCCT'
##single telo repeat looks like it's too short to use with just dumb matching.
##looking for repeat pairs, then multiply count later.

alltelos=tibble(chr=as.character(),
                length=as.numeric(),
                fwd=as.numeric(),
                rev=as.numeric(),
                asm=as.character())

for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes')
    prefixes=sub('\\.fasta$', '', list.files(genomedir, '.fasta$'))
    for (prefix in prefixes) {
        asmfile=file.path(genomedir, paste0(prefix, '.fasta'))
        telos=telocheck(asmfile, fwdtelo, revtelo) %>%
            mutate(fwd=fwd*2) %>%
            mutate(rev=rev*2) %>%
            mutate(asm=prefix)
        alltelos=bind_rows(alltelos, telos)
    }
}

fulltiginfo=inner_join(alltelos, alltiginfo, by=c('asm', 'length','chr'='tigname'))
tiginfofile=file.path(dbxdir, 'contig_info.csv')
write_csv(fulltiginfo, tiginfofile)

allteloplots=tibble(telopos=as.integer(),
                    percpos=as.numeric(),
                    telo=as.character(),
                    chr=as.character(),
                    name=as.character())

for (i in strains) {
    genomedir=file.path(datadir, i, 'genomes')
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

teloplotfile=file.path(dbxdir, 'telos.pdf')
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

