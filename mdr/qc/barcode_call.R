library(tidyverse)
library(BMS)
library(cowplot)
library(doParallel)
source('barcode_plot_functions.R')
cl=makeCluster(8)
registerDoParallel(cl, cores=8)
clusterCall(cl, function() library(tidyverse))

datadir='/mithril/Data/Nanopore/projects/methbin/barcode/qc'
dbxdir='~/Dropbox/timplab_data/mdr/barcode'

testfile=file.path(datadir, 'nebdcm_barcodes_filtered_10_motifs.txt')
ctrlfile=file.path(datadir, 'neb11_barcodes_filtered_10_motifs.txt')

#change this accordingly when barcoder bug is fixed
##bc_cols=c('readname', 'chrname', 'start', 'end', 'mapq', 'chrname2', 'GATC', 'CCWGG', 'ATCGAT', 'GTCGAC', 'CTCCAG', 'CTKVAG')
bc_cols=c('readname', 'chrname', 'start', 'end', 'mapq', 'GATC', 'GANTC', 'CCWGG', 'GCNGC')
bcinfo=read_tsv(testfile, col_names=bc_cols) %>%
    select(-start, -end, -mapq)
ctrlinfo=read_tsv(ctrlfile, col_names=bc_cols) %>%
    select(-start, -end, -mapq)


barcode='0010'
##realizing that you can't really reduce this to a binary classification
##well you can, but then weird things start to happen idk
thresholds=seq(0, 5, .1)
confusions=foreach(i=1:length(thresholds), .combine=rbind) %dopar% {
    barcodedreads=call_part_avg(bcinfo, thresholds[i])
    controlreads=call_part_avg(ctrlinfo, thresholds[i])
    
    tp=sum(barcodedreads$barcode==barcode)
    fp=sum(controlreads$barcode==barcode)
    
    tpr=tp/dim(barcodedreads)[1]
    fpr=fp/dim(controlreads)[1]
    
    conf=data.frame(thresh=thresholds[i], tpr=tpr, fpr=fpr)
    return(conf)
}
confusions=as_tibble(confusions) %>%
    mutate(area=tpr*(1-fpr))
##looking like 1.7 is the best threshold for dcm with min 10 motifs. continue evaluation this way


samps=c('neb11', 'neb19', 'neb17', 'neb15', 'nebdcm')
bcpops=tibble(pops=as.numeric(),
              barcode=as.character(),
              samp=as.character())
readlist=c()
for (i in samps) {
    file=file.path(datadir, paste0(i, '_barcodes_filtered_10_motifs.txt'))
    info=read_tsv(file, col_names=bc_cols) %>%
        select(-start, -end, -mapq)
    barcodedreads=call_part_avg(info, 1.7) %>%
        rowwise() %>%
        filter(!readname %in% readlist)
    if (i=='neb11') {
        readlist=c(readlist, barcodedreads$readname)
    }
    pops=tibble(pops=as.numeric(table(barcodedreads$barcode)),
                barcode=names(table(barcodedreads$barcode)),
                samp=i)
    bcpops=bind_rows(bcpops, pops)
}


popplots=list()
for (i in samps) {
    pops=bcpops %>%
        filter(samp==i)
    plot=plot_pops(pops, i)
    popplots[[i]]=plot
}


poppdf=file.path(dbxdir, 'barcode_pops.pdf')
pdf(poppdf, h=10, w=21)
print(plot_grid(popplots[[1]], popplots[[2]], popplots[[3]],  popplots[[4]], popplots[[5]]))
dev.off()
