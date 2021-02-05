library(tidyverse)


dbxdir='~/Dropbox/timplab_data/mdr5/abricate/'



##load distance matrix from clustalw and compute tree tree

library(ape)

distfile='/uru/Data/Nanopore/projects/mdr/all/clustalw_dists.csv'

dists=read_csv(distfile, skip=1, col_names=FALSE)
colnames(dists)=c('genes', dists$X1)
dists= dists %>%
    gather(key, value, -genes) %>%
    spread(genes, value) %>%
    select(-key)
tree=nj(as.dist(dists))
treefile='/uru/Data/Nanopore/projects/mdr/all/amrtree.nwk'
write.tree(tree, treefile)




##treeio to slice?
library(treeio)
treefile='/uru/Data/Nanopore/projects/mdr/all/amrtree_test.nwk'
##had to get rid of single quotes in order to get the thing to read
tree=read.newick(treefile)
tblamr=as_tibble(tree)

##try collapsing down one and two levels
collapse1=tblamr %>%
    group_by(parent) %>%
    summarise(seq=label[1])








##try using janky sorting
dists=read_csv(distfile, skip=1, col_names=FALSE)
colnames(dists)=c('genes', dists$X1)
dists= dists %>%
    gather(key, value, -genes) %>%
    spread(genes, value)

groupnum=0
groups=tibble(
    group=as.character(),
    seq=as.character())
while ( dim(dists)[1] > 0 ) {
    cols=colnames(dists)
    sorted_dists=dists %>%
        arrange(desc(get(cols[2])))
    diffs=abs(diff(sorted_dists %>% pull(get(cols[2]))))
    ##required cutoff is <95% similarity, and a diff of >5%
    minpos=sum(sorted_dists[,2]>95)
    maxdiffs=which(diffs>5)
    lastseq=maxdiffs[which(maxdiffs>minpos)[1]]

    ###define the group
    groupname=paste0('group', as.character(groupnum))

    seqs=sorted_dists$key[1:lastseq]
    group=tibble(group=groupname, seqs=seqs)
    groups=rbind(groups, group)
    groupnum=groupnum+1

    dists=dists %>%
        filter(!key %in% seqs) %>%
        select(
    
