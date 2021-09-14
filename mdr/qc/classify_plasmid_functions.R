library(tidyverse)
library(multidplyr)

cluster=new_cluster(6)
cluster_library(cluster, 'tidyverse')


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

get_classifications_from_umap <- function(embedplas, embedchr, numnear) {
    ##take each read in embedplas, return classifications of nearest reads in the umap
    ##numnear is how many nearest reads to return
    tigclass=embedplas %>%
        rowwise() %>%
        do(get_distances(., embedchr))

    return(tigclass)
}
