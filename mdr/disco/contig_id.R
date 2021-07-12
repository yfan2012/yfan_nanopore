library(tidyverse)

datadir='/mithril/Data/Nanopore/projects/methbin/disco'
dbdir='/atium/Data/ref/bacteria'

####look at blast results

##read accession key
faheaderfile=file.path(dbdir, 'all_bacteria_refs_faheaders.txt')
accinfo=read_table2(faheaderfile, col_names=c('acc', 'genus', 'species'))

##read blast results
blastfile=file.path(datadir, 'blast_contigs', 'flye_blast.tsv')
blastcols=c('tig', 'ref', 'ident', 'alen', 'mismatch', 'gaps', 'tigstart', 'tigend', 'refstart', 'refend', 'eval', 'bit')

blast=read_tsv(blastfile, col_names=blastcols, comment='#') %>%
    select(-c(mismatch, gaps, refstart, refend, eval, bit)) %>%
    filter(alen>1000 & ident>90)

##block=blast %>% filter(tig=='contig_1001' & ref=='NZ_WUQX01000001.1')
merge_overlaps <- function(block) {
    ##merge overlapping tig ranges
    ##return total bases covered on tig by the ref
    regions=tibble(start=block$tigstart,
                   end=block$tigend)
    numregions=dim(regions)[1]
    newnum=0
    
    while (newnum!=numregions) {
        newregions=regions[1,]
        if (dim(regions)[1]>=2) {
            for (i in 2:dim(regions)[1]) {
                newstart=regions$start[i]
                newend=regions$end[i]
                
                ##check for start or end expansion
                startexpand=which(newstart<newregions$start & newend>newregions$start)
                endexpand=which(newend>newregions$end & newstart<newregions$end)
                if (length(startexpand)>0) {
                    newregions$start[startexpand]=newstart
                }
                if (length(endexpand)>0) {
                    newregions$end[endexpand]=newend
                }

                ##check containment
                contain=which(newstart>=newregions$start & newend<=newregions$end)
                if (length(contain)<1) {
                    newregions=bind_rows(newregions, tibble(start=newstart, end=newend))
                }
            }
        }

        numregions=dim(regions)[1]
        newnum=dim(distinct(newregions))[1]
        regions=distinct(newregions)
    }
    return(regions)
}

merged=blast %>%
    group_by(tig, ref) %>%
    do(merge_overlaps(.)) %>%
    rowwise() %>%
    mutate(place=which(accinfo$acc==ref)) %>%
    mutate(g=accinfo$genus[place], s=accinfo$species[place])
 
mergecollapse=merged %>%
    mutate(numbases=end-start)%>%
    group_by(tig, g, s) %>%
    summarise(totalbases=sum(numbases))



####taxonomy stuff
##load taxids
taxidfile=file.path(dbdir, 'taxonomy', 'names_clean.dmp')
taxcols=c('id', 'name', 'specific', 'note', 'idk')
taxid=read_csv(taxidfile, col_names=taxcols) %>%
    filter(note=='scientific name') %>%
    ##apparently there is a genus of stick insect also called Bacillus. filter that out
    filter(id!=55087) %>%
    select(-c(idk, specific))
##load tax tree
treefile=file.path(dbdir, 'taxonomy', 'nodes_clean.dmp')
nodecols=c('tid', 'parent', 'rank', 'embl', 'divid', 'inheritdiv', 'genid', 'GC', 'mito', 'mgc', 'gb', 'hidden', 'comment', 'idk')
tree=read_csv(treefile, col_names=nodecols) %>%
    select(tid, parent)
    
idinfo=mergecollapse %>%
    rowwise() %>%
    mutate(taxid=case_when(length(which(g==taxid$name))==1 ~ taxid$id[which(g==taxid$name)][1],
                       length(which(g==taxid$name))==0 & length(which(paste0(g, ' ', s)==taxid$name))==1
                       ~ taxid$id[which(paste0(g, ' ', s)==taxid$name)][1]))

followback <- function(tid, tree) {
    ids=c(tid)
    lastid=tail(ids, n=1)
    newid=0
    while (newid!=lastid) {
        lastid=tail(ids, n=1)
        newid=tree$parent[which(tree$tid==lastid)]
        ids=c(ids, newid)
    }
    return(ids)
}

find_lca <- function(block, tree) {
    ##find lowest common ancestor of contig
    taxlist=followback(block$taxid[1], tree)
    if (dim(block)[1]>1) {
        for (i in 2:length(block$taxid)) {
            newlist=followback(block$taxid[i], tree)
            keepafter=which(taxlist %in% newlist)[1]
            end=length(taxlist)
            taxlist=taxlist[keepafter:end]
        }
    }
    result=tibble(tig=block$tig[1], id=taxlist[1])
    return(result)
}

classified=idinfo %>%
    group_by(tig) %>%
    do(find_lca(.,tree)) %>%
    rowwise() %>%
    mutate(name=taxid$name[which(taxid$id==id)])

classkeyfile=file.path(datadir, 'blast_contigs', 'classify_keys.csv')
write_csv(classified, classkeyfile)
                       
