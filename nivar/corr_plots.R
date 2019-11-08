library(ggplot2)
library(tidyverse)
library(tidyr)
library(gridExtra)
library(doParallel)
library(R.utils)
library(VennDiagram)
library(UpSetR)
library(ggpubr)

dbxdir='~/Dropbox/yfan/nivar/'
datadir='/uru/Data/Nanopore/projects/nivar/'

countsdir=paste0(dbxdir, 'motif_enrich/counts/')

##mismatch,ins,del per read
ciginfo=tibble(
    align.len = numeric(), 
    mismatch = numeric(),
    insert = numeric(),
    delete = numeric(),
    pore = character())
for (i in c('r9', 'r10')) {
    csv=paste0(datadir, 'align/', i, '/reference_' , i, '.md.sorted.csv')
    cig=read_csv(csv) %>%
        rowwise() %>%
        select(align.len, mismatch, insert, delete) %>%
        mutate(mismatch=mismatch/align.len) %>%
        mutate(insert=insert/align.len) %>%
        mutate(delete=delete/align.len) %>%
        mutate(pore=i)
    ciginfo=bind_rows(ciginfo, cig)
}
ciginfo=gather(ciginfo, key, value, -pore) %>%
    filter(key!='align.len')
pdf(paste0(dbxdir, 'read_errors.pdf'), width=12, height=8)
ggplot(ciginfo, aes(x=pore, y=value, colour=key, fill=key, alpha=.8)) +
    geom_violin() +
    xlab('pore') +
    ylim(0,.3) +
    theme_bw()
dev.off()



##compare number of corrections made
mumdir=paste0(datadir,'mummer/')
snpsfiles=list.files(mumdir, pattern='raw.snps', recursive=TRUE, include.dirs=TRUE)
rsnpsfiles=list.files(mumdir, pattern='raw_ref.snps', recursive=TRUE, include.dirs=TRUE)
snpsfiles=c(snpsfiles, rsnpsfiles)
corcounts=foreach(i=snpsfiles, .combine=rbind) %dopar% {
    numcors=countLines(paste0(mumdir, i))
    pore=strsplit(i, '_', fixed=TRUE)[[1]][1]
    corralg=strsplit(i, '_', fixed=TRUE)[[1]][2]
    if (corralg=='raw') {
        corralg='ref'
    }
    return(data.frame(samp=i, cors=numcors, pore=pore, alg=corralg))
}

totalplot=paste0(dbxdir, '/num_corrections.pdf')
pdf(totalplot, height=8.5, width=11)
ggplot(corcounts, aes(x=pore, y=cors)) +
    geom_bar(aes(fill=alg, colour=alg, alpha=.8),width=.5, stat='identity', position='dodge') +
    ggtitle('Corrections Made') +
    theme_bw()
dev.off()


##nearest neighbor error
distances=read_csv(paste0(dbxdir,'neighbor.csv'))
pdf(paste0(dbxdir,'neighbor.pdf'), height=8.5, width=11)
ggplot(distances, aes(x=distance)) +
    geom_histogram(data=distances[distances$pore=='r9',], fill='red', colour='red', alpha=.3, bins=40) +
    geom_histogram(data=distances[distances$pore=='r10',], fill='blue', colour='blue', alpha=.3, bins=40) +
    xlim(-10,100) +
    ggtitle('Nearest error in opposite pore') +
    xlab('Distance') +
    theme_bw()
dev.off()


##make tables of kmer populations
countfiles=list.files(countsdir, 'csv')
for (i in countfiles) {
    counts=read_csv(paste0(countsdir, i), col_names=c('seq', 'freq'))
    counts=counts %>%
        arrange(-freq) %>%
        mutate(rank=c(1:dim(counts)[1]))
pp
    name=substring(i, 1, nchar(i)-4)

    pdf(paste0(countsdir, name, '.pdf'), height=8.5, width=11)
    kmerhist=ggplot(counts, aes(x=freq)) +
        ggtitle('Error kmer histogram') +
        geom_histogram() +
            xlab('Number of errors') +
            theme_bw()
    rarefaction=ggplot(counts, aes(x=rank, y=freq)) +
        ggtitle('Error rarefaction') +
        geom_point() +
            xlab('Rank') +
            ylab('Frequency')+
            theme_bw()

    print(kmerhist)
    print(rarefaction)
    dev.off()
    
    pdf(paste0(countsdir, name, '_table.pdf'))
    print(grid.table(counts[1:20,1:2], rows=NULL))
    dev.off()
}


##plot a venn diagram of genomic locations that can be fixed
venn=read_csv(paste0(dbxdir, 'error_venn.csv'),skip=1, col_names=c('tig', 'both1', 'both2', 'r9', 'r10'))
venninfo=colSums(venn[2:5])
pdf(paste0(dbxdir, 'error_venn.pdf'), width=11, height=8.5)
venn.diagram(area1=venninfo[3]+venninfo[2],
             area2=venninfo[4]+venninfo[2],
             cross.area=venninfo[2],
             category=c('r9', 'r10'),
             lty=c('blank', 'blank'), 
             fill=c('light blue', 'pink'),
             alpha=c(.5, .5),
             scaled=TRUE)
dev.off()

for (i in c('r9', 'r10') ){
    vennr9=read_csv(paste0(dbxdir, i,'_correction_venn.csv'), skip=1, col_names=c('tig', 'all', 'fb_p', 'fb_r', 'r_p', 'fb', 'p', 'r'))
    vennr9info=unname(colSums(vennr9[2:8]))
    pdf(paste0(dbxdir, i, '_corr_ven.pdf'), width=11, height=8.5)
    expressionInput=c(fb=vennr9info[5],
                       pi=vennr9info[6],
                       ra=vennr9info[7],
                       `fb&pi`=vennr9info[2],
                       `fb&ra`=vennr9info[3],
                       `ra&pi`=vennr9info[4],
                      `fb&ra&pi`=vennr9info[1])
    print(upset(fromExpression(expressionInput)))
    dev.off()
}


##correction potential
corrs=c('pilon', 'racon', 'freebayes')

allranks=tibble(
    topmers = numeric(),
    common = numeric(),
    common_perc = numeric(),
    correction = character(), 
    r9perc = numeric(),
    r10perc = numeric()
)

corrpot=tibble(
    seq = character(),
    diffrank = numeric(),
    diffperc = numeric(),
    correction = character()
)
    
for (i in corrs) {
    r9=paste0(countsdir, '/nivar_r9_', i, '_raw.csv')
    r10=paste0(countsdir, '/nivar_r10_', i, '_raw.csv')

    r9counts=read_csv(r9, col_names=c('seq', 'freq')) %>%
        arrange(-freq) %>%
        mutate(perc_errors=freq/sum(freq)) %>%
        mutate(rank=1:length(freq)) %>%
        mutate(totalerr=1-cumsum(perc_errors))
    
    r10counts=read_csv(r10, col_names=c('seq', 'freq')) %>%
        arrange(-freq) %>%
        mutate(perc_errors=freq/sum(freq)) %>%
        mutate(rank=1:length(freq)) %>%
        mutate(totalerr=1-cumsum(perc_errors))

    pot9=r9counts %>%
        arrange(seq)
    
    pot10=r10counts %>%
        arrange(seq) %>%
        mutate(diffrank=abs(rank-pot9$rank)) %>%
        mutate(diffperc=abs(perc_errors-pot9$perc_errors)) %>%
        mutate(correction=i)

    corrpot=bind_rows(pot10[,c(1,6,7,8)])
    
    ranks=tibble(topmers=1:4096)
    ranks=ranks %>%
        rowwise() %>%
        mutate(common=sum(r9counts$seq[1:topmers] %in% r10counts$seq[1:topmers])) %>%
        mutate(common_perc=common/topmers) %>%
        mutate(correction=i) %>%
        add_column(r9perc=r9counts$totalerr[1:4096]) %>%
        add_column(r10perc=r10counts$totalerr[1:4096])
    
    allranks=bind_rows(allranks, ranks)
}    

pdf(paste0(dbxdir, 'motif_enrich/common_kmers.pdf'), height=8.5, width=11)
plot=ggplot(allranks, aes(x=topmers, y=common_perc, colour=correction)) +
    geom_line() +
    geom_line(linetype='dotted', aes(x=topmers, y=r9perc)) +
    geom_line(linetype='twodash',aes(x=topmers, y=r10perc)) +
    ggtitle('Erroneous 6-mers  in common between R9 and R10 assemblies') +
    xlab('6-mer rank') +
    ylab('Percent 6-mers in common') +
    ylim(0,1) +
    theme_bw()
print(plot)
dev.off()

pdf(paste0(dbxdir, 'motif_enrich/correction_potential.pdf'), height=12, width=8)
rankplot=ggplot(corrpot, aes(x=diffrank)) +
    geom_histogram(alpha=.8, colour='black') +
    ggtitle('Differences in 6-mer rank wrt error') +
    xlab('Difference in 6mer rank between r9 and r10') +
    ylab('Count') +
    theme_bw()
percplot=ggplot(corrpot, aes(x=diffperc)) +
    geom_histogram(bins=100, alpha=.8, colour='black') +
    ggtitle('Differences in how much error each 6-mer accounts for ') +
    xlab('Difference in percent error accounted for in r9 vs r10') +
    ylab('Count') +
    scale_y_log10() +
    xlim(0,.02) +
    theme_bw()
ggarrange(rankplot ,
          percplot,
          ncol=1, nrow=2, align='v')
dev.off()



##table classifying error kmers
kmerfeat=tibble(
    samp = character(), 
    gc = numeric(),
    hp = numeric(),
    entropy = numeric()
)

files=list.files(countsdir, '_ref.csv')
for (i in files) {
    path=paste0(dbxdir, 'motif_enrich/counts/', i)
    info=read_csv(path, col_names=c('seq', 'freq')) %>%
        rowwise() %>%
        mutate(hp=grepl('AAAA',seq) || grepl('TTTT',seq) || grepl('CCCC',seq) || grepl('GGGG',seq)) %>%
        mutate(gc=freq*(str_count(seq, 'G') + str_count(seq, 'C'))/6) %>%
        mutate(hpfreq=hp*freq) %>%
        mutate(plogp=freq/sum(freq)) %>%
        mutate(plogp=replace_na(plogp,0))

    plogp=log2(info$freq/sum(info$freq))*info$freq/sum(info$freq)
    plogp[is.na(plogp)] = 0
    entropy=-sum(plogp)
    gc=sum(info$gc)/sum(info$freq)
    hp=sum(info$hpfreq)/sum(info$freq)

    kmerfeat=bind_rows(kmerfeat, tibble(samp=i, gc=gc, hp=hp, entropy=entropy))
}
kmerfeat=kmerfeat %>%
    rowwise() %>%
    mutate(samp=paste0(strsplit(samp, '_', fixed=TRUE)[[1]][2], '_', strsplit(samp, '_', fixed=TRUE)[[1]][3])) %>%
    mutate(pore=strsplit(samp, '_', fixed=TRUE)[[1]][1]) %>%
    mutate(alg=strsplit(samp, '_', fixed=TRUE)[[1]][2])

pdf(paste0(dbxdir, 'kmerfeat.pdf'), height=12, width=8)
gc=ggplot(kmerfeat, aes(x=pore, y=gc)) +
    geom_bar(aes(fill=alg, colour=alg, alpha=.8), width=.5, stat='identity', position='dodge') +
    xlab('Pore') +
    ylab('% GC content of error sequences') +
    ggtitle('GC Content') +
    theme_bw()
hp=ggplot(kmerfeat, aes(x=pore, y=hp)) +
    geom_bar(aes(fill=alg, colour=alg, alpha=.8), width=.5, stat='identity', position='dodge') +
    xlab('Pore') +
    ylab('% homopolymers of error sequences') +
    ggtitle('Homopolymer Content') +
    theme_bw()
ent=ggplot(kmerfeat, aes(x=pore, y=entropy)) +
    geom_bar(aes(fill=alg, colour=alg, alpha=.8), width=.5, stat='identity', position='dodge') +
    xlab('Pore') +
    ylab('Entropy of error sequences') +
    ggtitle('Entropy') +
    theme_bw()
ggarrange(gc + rremove('x.text') + rremove('x.title'),
          hp + rremove('x.text') + rremove('x.title'),
          ent,
          ncol=1, nrow=3, align='v')
dev.off()
pdf(paste0(dbxdir, 'kmerfeat_table.pdf'))
grid.table(kmerfeat[1:8,])
dev.off()



##plot types of errors corrected
types=read_csv(paste0(dbxdir, 'correction_types.csv')) %>%
    gather(key, value, -pore, -alg)
pdf(paste0(dbxdir, 'correction_types.pdf'), height=12, width=8)
fb=ggplot(types[types$alg=='freebayes',], aes(x=pore, y=value)) +
    geom_bar(aes(fill=key, colour=key, alpha=.8), width=.5, stat='identity', position='dodge') +
    xlab('Pore') +
    ggtitle('Freebayes') +
    theme_bw()
pilon=ggplot(types[types$alg=='pilon',], aes(x=pore, y=value)) +
    geom_bar(aes(fill=key, colour=key, alpha=.8), width=.5, stat='identity', position='dodge') +
    xlab('Pore') +
    ggtitle('Pilon') +
    theme_bw()
racon=ggplot(types[types$alg=='racon',], aes(x=pore, y=value)) +
    geom_bar(aes(fill=key, colour=key, alpha=.8), width=.5, stat='identity', position='dodge') +
    xlab('Pore') +
    ggtitle('Racon') +
    theme_bw()
allcors=ggplot(corcounts, aes(x=pore, y=cors)) +
    geom_bar(aes(fill=alg, colour=alg, alpha=.8),width=.5, stat='identity', position='dodge') +
    ggtitle('Corrections Made') +
    theme_bw()
ggarrange(allcors + rremove('x.text') + rremove('x.title'),
          fb + rremove('x.text') + rremove('x.title'),
          pilon+ rremove('x.text') + rremove( 'x.title'),
          racon,  
          ncol=1, nrow=4, align='v')
dev.off()

    


##assembly stats and yield tables
asmfile=paste0(dbxdir, 'asm_stats.csv')
asm=read_csv(asmfile,
             col_names=c('asm', 'num','n50', 'longest', 'shortest', 'total')) %>%
    rowwise() %>%
    mutate(asm=strsplit(asm, '/', fixed=TRUE)[[1]][8])
pdf(paste0(dbxdir, 'asm_stats.pdf'))
grid.table(asm)
dev.off()

yieldfile=paste0(dbxdir, 'yields.csv')
yield=read_csv(yieldfile,
               col_names=c('file', 'yield', 'numreads', 'avglen')) %>%
    rowwise() %>%
    mutate(file=strsplit(file, '/', fixed=TRUE)[[1]][8]) %>%
    mutate(yield=round(yield/1000000, 2)) %>%
    mutate(numreads=round(numreads/1000000, 2)) %>%
    group_by(file) %>%
    summarise_each(funs(sum, mean), yield, numreads, avglen)
yield=yield[c(1,2,3,7)]
colnames(yield)=c('file','yield','numreads','avglen')

pdf(paste0(dbxdir, 'yields.pdf'))
grid.table(yield)
dev.off()
