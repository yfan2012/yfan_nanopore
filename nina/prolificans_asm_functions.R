library(tidyverse)
library(Biostrings)
library(RColorBrewer)

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
    
    teloinfo=tibble(chr=names(asm)) %>%
        mutate(length=width(asm[chr])) %>%
        rowwise() %>%
        mutate(fwdlow=sum(str_locate_all(as.character(asm[chr]), fwdtelo)[[1]][,1] < 1000)) %>%
        mutate(fwdhigh=sum(str_locate_all(as.character(asm[chr]), fwdtelo)[[1]][,1] > length-1000)) %>%
        mutate(revlow=sum(str_locate_all(as.character(asm[chr]), revtelo)[[1]][,1] < 1000)) %>%
        mutate(revhigh=sum(str_locate_all(as.character(asm[chr]), revtelo)[[1]][,1] > length-1000)) %>%
        mutate(fwd=fwdlow+fwdhigh) %>%
        mutate(rev=revlow+revhigh)
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

library(BSgenome)
library(GenomicRanges)

suggest_mito_breaks <- function(coordsfile, asmtiginfo, asmfile, outfile) {
    ##find best mito candidate and write out to file
    ##best candidate: longest alignment to a tig with no telos
    coordcols=c('rstart', 'rend', 'qstart', 'qend', 'raln', 'qaln','ident','rlen', 'qlen', 'rfrac', 'qfrac', 'ref', 'qry')
    notelos=asmtiginfo %>%
        filter(fwd==0 & rev==0)
    coords=read_tsv(coordsfile, col_names=coordcols, col_types=cols(), progress=FALSE) %>%
        rowwise() %>%
        filter(ref %in% notelos$chr)
    asm=readDNAStringSet(asmfile)
    best=coords[coords$qfrac==max(coords$qfrac),]

    mitotigs=c()
    if (dim(best)[1]>1) {
        ##if more than one best align, check for same chr
        mitotigs=unique(best$ref)
        if (length(mitotigs)>1) {
            print('more than one mito tig? taking best and deleting the rest')
        }
    }
    best=best[best$ident==max(best$ident)[1],]
        
    mitoregion=GRanges(seqnames=best$ref, ranges=IRanges(start=best$rstart, end=best$rend))
    mitoseq=getSeq(asm, mitoregion)
    names(mitoseq)=paste0(notelos$asm[1], '_mito')
    newasm=c(asm[names(asm)[names(asm)!=best$ref & !(names(asm) %in% mitotigs)]], mitoseq)
    writeXStringSet(newasm, outfile, format='fasta')
}

    
