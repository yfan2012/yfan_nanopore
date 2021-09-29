library(tidyverse)
library(multidplyr)

#cluster=new_cluster(6)
#cluster_library(cluster, 'tidyverse')


get_distances <- function(plasread, embedchr) {
    colnames=names(table(embedchr$label))
    classcounts=as_tibble(matrix(0, nrow=2, ncol=length(colnames)))[1,]
    names(classcounts)=colnames
    
    distances=embedchr %>%
        mutate(xdist=(plasread$x-x)^2) %>%
        mutate(ydist=(plasread$y-y)^2) %>%
        mutate(euclidian=sqrt(xdist+ydist)) %>%
        arrange(euclidian)
    top=distances[1:50,] %>%
        group_by(label) %>%
        summarise(counts=n())
    for (i in top$label) {
        classcounts[1,names(classcounts)==i]=top$counts[top$label==i]
    }
    plasclass=bind_cols(plasread, classcounts)
    return(plasclass)
}

classify_umap_neighbors <- function(embedplas, embedchr, numnear) {
    ##take each read in embedplas, return classifications of nearest reads in the umap
    ##numnear is how many nearest reads to return
    tigclass=embedplas %>%
        rowwise() %>%
        do(get_distances(., embedchr))

    return(tigclass)
}


classify_plasmid_reads_distance <- function(plasinforow, tigsbulk) {
    ##take a plasmid read
    ##take bulk info
    ##get euclidean distance from each tig
    ##make sure everything is already scaled
    plasinforow=as_tibble(plasinforow)
    colnames=c(tigsbulk$chrname)
    dists=matrix(0, 1, length(colnames))
    colnames(dists)=colnames

    for (i in 1:dim(tigsbulk)[1]) {
        tigsbulkrow=tigsbulk[i,]
        rowdata=plasinforow %>%
            select(-readname, -chrname)
        tigsbulkdata=tigsbulkrow %>%
            select(-chrname)
        pair=bind_rows(rowdata, tigsbulkdata)
        distance=dist(pair)
        
        dists[1,tigsbulkrow$chrname]=distance
    }

    distinfo=tibble(chrname=plasinforow$chrname)
    distinfo=bind_cols(distinfo, as_tibble(dists))
    
    return(distinfo)
}
