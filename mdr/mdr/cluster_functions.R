filter_reads <- function(fullbcinfo, fullbccounts, maxna, mincount, maxreads, chrkey, cluster) {
    cluster_copy(cluster, 'maxreads')
    
    nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
    lowna=nacount[nacount<maxna]
    keepmotifs=names(lowna)
    
    bcinfo=fullbcinfo %>%
        select(all_of(keepmotifs))
    
    
    ##this set of motifs is really robust to high motif count requirement apparently
    ##may need to knowck this back down in order to relax the nacount requirement later
    bccounts=fullbccounts %>%
        select(all_of(keepmotifs)) %>%
        filter(across(c(-readname, -chrname), ~ .x>=mincount))
    cluster_copy(cluster, 'bccounts')

    
    ##filter out reads with insufficient motif coverage
    bcfilt=bcinfo %>%
        filter(readname %in% bccounts$readname) %>%
        mutate(chrname=case_when(substr(chrname,1,7)=='NZ_JXXK' ~ 'NZ_JXXK',
                                 substr(chrname,1,7)=='NZ_KE13' ~ 'NZ_KE13',
                                 TRUE ~ chrname)) %>%
        rowwise() %>%
        mutate(species=chrkey$species[chrkey$acc==chrname]) %>%
        mutate(type=chrkey$type[chrkey$acc==chrname]) %>%
        group_by(chrname) %>%
        partition(cluster)

    
    bcfiltsub=bcfilt %>%
        do(checkfilt(., maxreads, bccounts)) %>%
        collect() %>%
        ungroup() %>%
        filter(complete.cases(.))
    bcdata=bcfiltsub %>%
        select(-chrname, -readname, -species, -type)
    
    bcumap=umap(bcdata)
    embed=tibble(x=bcumap$layout[,1],
                 y=bcumap$layout[,2],
                 label=bcfiltsub$species,
                 type=bcfiltsub$type,
                 read=bcfiltsub$readname)
    return(embed)
}

plotembed <- function(embed, clusterfile, alpha) {
    myshapes=c(3,2,1)

    embedchr=embed %>%
        filter(type=='chr')
    embedplas=embed %>%
        filter(type=='plas')
    
    pdf(clusterfile, h=9, w=13)
    plot=ggplot(embedchr, aes(x=x, y=y, colour=label))+
        geom_point(alpha=alpha, size=.5) +
        scale_colour_brewer(palette='Set2') +
        theme_bw()
    mainplot=plot +
        geom_point(data=embedplas, aes(x=x, y=y, shape=label), inherit.aes=FALSE) +
        scale_shape_manual(values=myshapes) +
        theme(legend.position = 'none')
    print(mainplot)
    sep=plot +
        facet_wrap(~label) +
        geom_point(data=embedplas, aes(x=x, y=y, shape=label), size=.5, alpha=alpha, inherit.aes=FALSE) +
        scale_shape_manual(values=myshapes) +
        theme(legend.position = "none")
    print(sep)
    dev.off()
}



filter_align_reads <- function(fullbcinfo, fullbccounts, paffilt, maxna, mincount, maxreads, chrkey, cluster) {
    cluster_copy(cluster, 'maxreads')
    
    nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
    lowna=nacount[nacount<maxna]
    keepmotifs=names(lowna)
    
    bcinfo=fullbcinfo %>%
        select(all_of(keepmotifs)) %>%
        filter(readname %in% paffilt$qname)
    
    
    ##this set of motifs is really robust to high motif count requirement apparently
    ##may need to knowck this back down in order to relax the nacount requirement later
    bccounts=fullbccounts %>%
        select(all_of(keepmotifs)) %>%
        filter(across(c(-readname, -chrname), ~ .x>=mincount))
    cluster_copy(cluster, 'bccounts')

    
    ##filter out reads with insufficient motif coverage
    bcfilt=bcinfo %>%
        filter(readname %in% bccounts$readname) %>%
        mutate(chrname=case_when(substr(chrname,1,7)=='NZ_JXXK' ~ 'NZ_JXXK',
                                 substr(chrname,1,7)=='NZ_KE13' ~ 'NZ_KE13',
                                 TRUE ~ chrname)) %>%
        rowwise() %>%
        mutate(species=chrkey$species[chrkey$acc==chrname]) %>%
        mutate(type=chrkey$type[chrkey$acc==chrname]) %>%
        group_by(chrname) %>%
        partition(cluster)

    
    bcfiltsub=bcfilt %>%
        do(checkfilt(., maxreads, bccounts)) %>%
        collect() %>%
        ungroup() %>%
        filter(complete.cases(.))
    bcdata=bcfiltsub %>%
        select(-chrname, -readname, -species, -type)
    
    bcumap=umap(bcdata)
    embed=tibble(x=bcumap$layout[,1],
                 y=bcumap$layout[,2],
                 label=bcfiltsub$species,
                 type=bcfiltsub$type,
                 read=bcfiltsub$readname)
    return(embed)
}
