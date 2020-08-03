library(tidyverse)
library(UpSetR)


dbxdir='~/Dropbox/timplab_data/mdr5/abricate/'


##load combined abricate report and make an upset plot
mergefile='/uru/Data/Nanopore/projects/mdr/all/abricate_all.tsv'

amr=read_tsv(mergefile, col_types = cols(.default = "c")) %>%
    select(-NUM_FOUND) 
##changes to 0 and 1
amr[, -1][amr[, -1]=='.']='0'
amr[, -1][amr[, -1]!='0']='1'
amr[, -1]=sapply(amr[, -1], as.integer)


amr=amr %>%
    mutate(file=gsub('/uru/Data/Nanopore/projects/mdr/', '', `#FILE`)) %>%
    mutate(file=gsub('.all.tsv', '', file)) %>%
    select(-`#FILE`) 

amr=amr %>%
    gather(variable, value, -file) %>%
    spread(file, value)

outfile=paste0(dbxdir, 'amr_upset.pdf')
pdf(outfile, h=7, w=15)
sets=colnames(amr[-1])
upset(as.data.frame(amr), sets=sets, order.by = "freq", empty.intersections = "on")
dev.off()


##just some exploration of what's diffent between the native samp and the other two
native_extra=amr %>%
    filter(`illumina/amr/MDRstool_16_metaspades`==0) %>%
    filter(`MDRstool_16/amr/MDRstool_16_pcr_100m`==0)
native_missing=amr %>%
    filter(`illumina/amr/MDRstool_16_metaspades`==1) %>%
    filter(`MDRstool_16/amr/MDRstool_16_pcr_100m`==1) %>%
    filter(`MDRstool_16/amr/MDRstool_16_native_100m`==0)



##load originals
datadir='/uru/Data/Nanopore/projects/mdr/MDRstool_16/amr/'
files=list.files(datadir, patter='all.tsv')

library(GenomicRanges)

grouphits <- function(abrifile, prefix){
    ###group hits by contig locations
    colnames=c('file','tig','start','end','strand','gene','cov','covmap','gaps','perccov','ident','db','acc','product','resistance')
    abri=read_tsv(abrifile, col_names=colnames)

    groups=tibble(
        group=as.character(),
        accs=as.character())

    gr=GRanges(seqnames=abri$tig,
               ranges=IRanges(start=abri$start, end=abri$end, names=abri$acc))


    for (i in 1:length(gr)) {
        q=gr[i]
        overlaps=findOverlaps(q,gr)
        ##should at least hit to itself, so shouldn't need to account for no indexes
        indexes=subjectHits(overlaps)
        groupacc=paste(sort(names(gr[indexes])), collapse=',')

        groupname=paste0(prefix, i)
        groups=rbind(groups, tibble(group=groupname, accs=groupacc))
    }
            
    return(groups)
}

illfile=paste0(datadir, files[1])
illgroups=grouphits(illfile, 'ill')
natfile=paste0(datadir, files[2])
natgroups=grouphits(natfile, 'nat')
pcrfile=paste0(datadir, files[3])
pcrgroups=grouphits(pcrfile, 'pcr')

allgroups=rbind(illgroups, natgroups, pcrgroups)

groupinfo=tibble(
    group=as.character(),
    ill=as.integer(),
    nat=as.integer(),
    pcr=as.integer())

##this does not feel elegant, but i am suffering from menstrual pain atm, so will just go with it for now. 
uniquegroups=unique(allgroups$accs)
for (i in uniquegroups) {
    ##check what it's in
    ill=0
    nat=0
    pcr=0
    if (i %in% illgroups$accs) {
        ill=1
    }
    if (i %in% natgroups$accs) {
        nat=1
    }
    if (i %in% pcrgroups$accs) {
        pcr=1
    }
    newinfo=tibble(group=i, ill=ill, nat=nat, pcr=pcr)
    groupinfo=rbind(groupinfo, newinfo)
}

outfile=paste0(dbxdir, 'amrgroups_upset.pdf')
pdf(outfile, h=7, w=15)
sets=colnames(groupinfo[-1])
upset(as.data.frame(groupinfo), sets=sets, order.by = "freq", empty.intersections = "on")
dev.off()


##investigate some outgroups
natout=groupinfo %>%
    filter(ill==0 & nat==1 & pcr==0)
illout=groupinfo %>%
    filter(ill==1 & nat==0 & pcr==0)
pcrout=groupinfo %>%
    filter(ill==0 & nat==0 & pcr==1)
