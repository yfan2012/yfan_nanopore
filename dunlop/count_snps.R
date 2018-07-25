library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)

vcf='~/Dropbox/yfan/dunlop/parsnp/pilon17_parsnp/pilon17.vcf'

snptab=read_tsv(vcf, comment='##') %>%
    filter(FILTER=='PASS')


snps=matrix(nrow=4, ncol=4)
for (i in 1:4) {
    for (j in 1:4) {
        snps[i, j]=sum(snptab[,i+9]!=snptab[,j+9])
    }
}

names=colnames(snptab)[10:13]


for (i in 1:length(names)){
    newname=gsub('.fasta', '',names[i])
    names[i]=newname
}

rownames(snps)=names
colnames(snps)=names

pdf('/home/yfan/Dropbox/yfan/dunlop/numsnps.pdf', width=10, height=5)
grid.table(snps)
dev.off()

