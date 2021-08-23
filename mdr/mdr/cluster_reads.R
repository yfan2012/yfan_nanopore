library(tidyverse)
library(umap)
library(multidplyr)
source('~/Code/yfan_nanopore/mdr/qc/barcode_plot_functions.R')
source('~/Code/yfan_nanopore/mdr/mdr/cluster_functions.R')

cluster=new_cluster(6)
cluster_library(cluster, 'tidyverse')
cluster_copy(cluster, 'checkfilt')

datadir='/mithril/Data/Nanopore/projects/methbin/mdr'
dbxdir='~/gdrive/mdr/mdr'
prefix='200708_mdr_stool16native'

barcodelist=file.path('~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt')
motifinfo=read_tsv(barcodelist, col_names=FALSE)
bc_cols=c('readname', 'chrname', motifinfo$X1)

barcodedatafile=file.path(datadir, 'barcode', prefix, paste0(prefix, '_barcodes.txt'))
fullbcinfo=read_tsv(barcodedatafile, col_names=bc_cols, na=c('None'))
barcodecountsfile=file.path(datadir, 'barcode', prefix,paste0(prefix, '_barcodes_motifcounts.txt'))
fullbccounts=read_tsv(barcodecountsfile, col_names=bc_cols)

##chr vs plasmid labels
chrkeyfile=file.path(datadir, 'ref/mdr_refs.txt')
chrkey=read_csv(chrkeyfile)




####try .3 maxna, 5 mincounts, 10k maxreads
embed=filter_reads(fullbcinfo, fullbccounts, .3, 5, 10000, chrkey, cluster)
clusterfile=file.path(dbxdir, 'cluster_10k.pdf')
plotembed(embed, clusterfile, alpha=.2)

####try .3 maxna, 5 mincounts, 1k maxreads
embed=filter_reads(fullbcinfo, fullbccounts, .3, 5, 1000, chrkey, cluster)
clusterfile=file.path(dbxdir, 'cluster_1k.pdf')
plotembed(embed, clusterfile, alpha=.2)

####try .3 maxna, 5 mincounts, 500 maxreads
embed=filter_reads(fullbcinfo, fullbccounts, .3, 5, 500, chrkey, cluster)
clusterfile=file.path(dbxdir, 'cluster_500.pdf')
plotembed(embed, clusterfile, alpha=.2)





####look at individual clusters
embed=filter_reads(fullbcinfo, fullbccounts, .3, 5, 1000, chrkey, cluster)

##grab relevant reads
cluster1=embed %>%
    filter(y < -2)

nacount=colSums(is.na(fullbcinfo)/dim(fullbcinfo)[1])
lowna=names(nacount[nacount<.3])

clustfilt=fullbcinfo %>%
    filter(readname %in% cluster1$read) %>%
    select(all_of(lowna))

clustdata=clustfilt %>%
    select(-readname, -chrname)
clustumap=umap(clustdata)
clustembed=tibble(x=clustumap$layout[,1],
                  y=clustumap$layout[,2],
                  chrname=clustfilt$chrname) %>%
                  mutate(chrname=case_when(substr(chrname,1,7)=='NZ_JXXK' ~ 'NZ_JXXK',
                                           substr(chrname,1,7)=='NZ_KE13' ~ 'NZ_KE13',
                                           TRUE ~ chrname)) %>%
                  rowwise() %>%
                  mutate(label=chrkey$species[chrkey$acc==chrname]) %>%
                  mutate(type=chrkey$type[chrkey$acc==chrname])
clusterfile=file.path(dbxdir, 'subcluster1_1k.pdf')
plotembed(clustembed, clusterfile, alpha=.5)
                  





####pull out plasmids
embed=filter_reads(fullbcinfo, fullbccounts, .3, 5, 10000, chrkey, cluster)
klebplas=embed %>%
    filter(y<0) %>%
    filter(label=='kleb')
thetaplas=embed %>%
    filter(y>2) %>%
    filter(label=='theta')
plas=bind_rows(klebplas, thetaplas)

##read align
paffile=file.path(datadir, 'align', paste0(prefix, '.paf'))
pafcols=c('qname', 'qlen', 'qstart', 'qend', 'strand', 'rname', 'rlen', 'rstart', 'rend', 'match', 'alen', 'mapq', 'tp', 'cm', 's1', 's2', 'dv', 'rl')
paf=read_tsv(paffile, col_names=pafcols) %>%
    select(qname, qlen, qstart, qend, rname, rlen, rstart, rend, match, alen, mapq) %>%
    filter(qname %in% plas$read) %>%
    rename(read = qname)

plasinfo=full_join(paf, plas, by='read')
plascsv=file.path(dbxdir, 'cross_reads.csv')
write_csv(plasinfo, plascsv)







####filter align to use for read filter later
paffilt=read_tsv(paffile, col_names=pafcols) %>%
    group_by(qname) %>%
    filter(n()==1) %>%
    filter(mapq>59)



embed=filter_align_reads(fullbcinfo, fullbccounts, paffilt, .3, 5, 10000, chrkey, cluster)
clusterfile=file.path(dbxdir, 'cluster_10k_alignfilt.pdf')
plotembed(embed, clusterfile, alpha=.2)


theta=embed %>%
    filter(label=='theta') %>%
    filter(type=='plas')
thetacoords=paffilt %>%
    rowwise() %>%
    filter(qname %in% theta$read) %>%
    select(-tp, -cm, -s1, -s2, -dv, -rl) %>%
    rename(read=qname)
thetainfo=full_join(theta, thetacoords)

kleb=embed %>%
    filter(label=='kleb') %>%
    filter(type=='plas')
klebcoords=paffilt %>%
    rowwise() %>%
    filter(qname %in% kleb$read) %>%
    select(-tp, -cm, -s1, -s2, -dv, -rl) %>%
    rename(read=qname)
klebinfo=full_join(kleb, klebcoords)
