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




##load distance matrix from clustalw and plot tree
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





##try using janky sorting
dists=read_csv(distfile, skip=1, col_names=FALSE)
colnames(dists)=c('genes', dists$X1)
dists= dists %>%
    gather(key, value, -genes) %>%
    spread(genes, value)

gap=as.numeric()
while ( dim(dists)[1] > 0 ) {
    cols=colnames(dists)
    sorted_dists=dists %>%
        arrange(desc(get(cols[2])))
    diffs=abs(diff(sorted_dists %>% pull(get(cols[2]))))
    ##required cutoff is <95% similarity, and a drop of > 5%
