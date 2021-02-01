library(tidyverse)

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

cov=get_asm_cov(npcovfile, illcovfile)

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

tiginfo=get_tigs_info(cov)

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

tiginfofile=file.path(dbxdir, 'contig_info.csv')
write_csv(alltiginfo, tiginfofile)

breakinfofile=file.path(dbxdir, 'zero_cov.csv')
write_csv(allbreaks,breakinfofile)




            
