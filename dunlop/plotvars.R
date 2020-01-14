library(ggplot2)
library(tidyverse)

datadir='/uru/Data/dunlop/'
prefix='190918_dunlop'
dbxdir=paste0('~/Dropbox/yfan/dunlop/',prefix,'/')

varsdir=paste0(datadir,prefix,'/vars')
varfiles=list.files(varsdir, 'csv')

##read in all the vars data
vars=read_csv(paste0(varsdir, '/', varfiles[1]), col_names=c('chr', 'pos','ref', 'alt', 'numref', 'numalt', 'gene', 'prod', 'start', 'end')) %>%
    mutate(prefix=basename(varfiles[1]))
for (i in varfiles[-1]) {
    tmpvars=read_csv(paste0(varsdir, '/', i), col_names=c('chr', 'pos','ref', 'alt', 'numref', 'numalt', 'gene', 'prod', 'start', 'end')) %>%
        mutate(prefix=basename(i))
    vars=rbind(vars, tmpvars)
}


vars=vars %>%
    filter(numref>0) %>%
    group_by(prefix, chr, pos) %>%
    mutate(altperc=numalt/(numref+sum(numalt))) %>%
    mutate(refperc=numref/(numref+sum(numalt)))


##report all vars with more than 3 reads supporting it
highev=vars[vars$numalt>3,] %>%
    arrange(desc(pos))
    ##filter(gene!='unk')
write.table(data.frame(highev), paste0(dbxdir, 'filtered_muts.csv'), sep=',')

##info by position - this includes the unks
meanpos=vars %>%
    group_by(prefix, pos) %>%
    summarise(meanperc=sum)

##info by genes
geneinfopos=vars %>%
    filter(gene!='unk') %>%
    mutate(genelen=abs(as.numeric(end)-as.numeric(start))) %>%
    group_by(prefix, gene, pos) %>%
    summarise(cov=sum(numalt)+numref[1], numpos=length(pos), genelen=genelen[1]) %>%
    group_by(prefix, gene) %>%
    ##find number of unique positions that are called as a var by one (1) read
    summarise(altpos=length(unique(pos))/mean(cov)/genelen[1])
    ##summarise(altreads=sum(numalt)/mean(genelen))
geneinfomut=vars %>%
    filter(gene!='unk') %>%
    mutate(genelen=abs(as.numeric(end)-as.numeric(start))) %>%
    group_by(prefix, gene) %>%
    ##find number of unique positions that are called as a var by one (1) read
    ##summarise(altspos=length(unique(pos))/mean(genelen)) %>%
    summarise(altreads=sum(altperc)/mean(genelen))


genefile=paste0(dbxdir, 'genepos_dist.pdf')
pdf(genefile, height=8, width=15)
ggplot(geneinfopos, aes(x=prefix, y=altpos, fill=prefix, colour=prefix, alpha=.5)) +
    geom_violin(scale='width') +
    ggtitle('Mutation Loci per Gene (supported by 1 read, normalized by gene length and coverage)') +
    theme_bw()
ggplot(geneinfomut, aes(x=prefix, y=altreads, fill=prefix, colour=prefix, alpha=.5)) +
    geom_violin(scale='width') +
    ggtitle('Mutations per Gene (supported by 1 read, normalized by gene length and coverage)') +
    theme_bw()
dev.off()
   

              
