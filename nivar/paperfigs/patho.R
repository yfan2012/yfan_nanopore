library(tidyverse)
library(tidyr)
library(gridExtra)
library(RColorBrewer)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(R.utils)
library(foreach)
library(doParallel)
library(colorspace)

dbxdir='~/Dropbox/yfan/nivar/'
datadir='/uru/Data/Nanopore/projects/nivar/paperfigs/'

asmfile='/uru/Data/Nanopore/projects/nivar/paperfigs/assembly_final/nivar.final.fasta'
asmgff='/uru/Data/Nanopore/projects/nivar/paperfigs/annotation_final/nivar.final.gff'

reffile='/uru/Data/Nanopore/projects/nivar/reference/candida_nivariensis.fa'
glafile='/uru/Data/Nanopore/projects/nivar/reference/medusa_fungi/candida_glabrata.fa'
xufile='/uru/Data/Nanopore/projects/nivar/paperfigs/patho/glabrata_xu.fa'

regionsfile='/uru/Data/Nanopore/projects/nivar/paperfigs/patho/gpicwp_regions.csv'

##colors
colors=tibble(scale=seq(0, .8, .1)) %>%
    mutate(light=lighten('#8DA0CB', scale)) %>%
    mutate(dark=darken('#8DA0CB', scale))
colvec=c(rev(colors$light), colors$dark[-1])


##read in regions file from xu (asterisks removed)
regiontable=read_csv(regionsfile) %>%
    filter(!str_detect(Status, 'Removed')) %>%
    select(-Source, -Note, -CommonName)
regions=GRanges(seqnames=regiontable$Chrom,
                ranges=IRanges(start=regiontable$Start, end=regiontable$End),
                strand=regiontable$Direction,
                gene=regiontable$Systematic_Name)

##rename so the chrs match the table
xugla=readDNAStringSet(xufile)
chrnames=tibble(chr='Chr', letter=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M')) %>%
    mutate(chrname=paste0(chr, letter))
names(xugla)=chrnames$chrname

##get gene seqs
geneseqs=getSeq(xugla, regions)
names(geneseqs)=regiontable$Systematic_Name

##write to fasta
geneseqsfa=file.path(datadir, 'patho', 'gpicwp.fa')
writeXStringSet(geneseqs, geneseqsfa, format='fasta')





##compare gpicwp hits between ref and asm
cnames=c('gene', 'chr', 'ident', 'alignlen', 'mismatches', 'gap', 'qstart', 'qend', 'sstart', 'send', 'eval', 'store')
nivtsv=file.path(datadir, 'patho', 'nivar.final.gpicwp_hits.tsv')
reftsv=file.path(datadir, 'patho', 'candida_nivariensis.tsv')
nivhits=read_tsv(nivtsv, comment='#', col_names=cnames) %>%
    mutate(genome='asm')
refhits=read_tsv(reftsv, comment='#', col_names=cnames) %>%
    mutate(genome='ref')

genecounts <- function(hits) {
    ##count genes by weight
    minlen=width(geneseqs[hits$gene[1]])

    hits=hits %>%
        mutate(contribute=alignlen/minlen)

    refs=sum(hits$contribute[hits$genome=='ref'])
    asms=sum(hits$contribute[hits$genome=='asm'])
        
    genecount=tibble(gene=hits$gene[1],
                     asm=c('ref', 'asm'),
                     counts=c(refs, asms),
                     numhits=c(sum(hits$genome=='ref'), sum(hits$genome=='asm')))
    return(genecount)
}

allhits=rbind(nivhits, refhits) %>%
    group_by(gene) %>%
    do(genecounts(.)) %>%
    na.omit %>%
    rowwise() %>%
    mutate(status=regiontable$Status[regiontable$Systematic_Name==gene])
allhitscsv=file.path(dbxdir, 'paperfigs','raw', 'copynum.csv')
write_csv(allhits, allhitscsv)

allhitscorr=spread(allhits[,1:3], asm, counts)
allhitscorrcsv=file.path(dbxdir, 'paperfigs','raw', 'copynum_corr.csv')
write_csv(allhitscorr, allhitscorrcsv)

copynumpdf=file.path(dbxdir, 'paperfigs', 'raw', 'copynum.pdf')
pdf(copynumpdf, height=8, width=7)
plot=ggplot(allhits, aes(x=numhits, fill=asm, colour=asm, alpha=.3)) +
    geom_histogram(stat='count', binwidth=1) +
    facet_wrap(. ~ asm, ncol=1) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    xlim(0,10) +
    scale_x_continuous(name='Number of Hits', limits=c(0,10), breaks=seq(0,10,1)) +
    ggtitle('Telomeric Gene Copy Number') +
    theme_bw()
print(plot)
asmheatmap=ggplot(allhits, aes(numhits, counts)) +
    geom_bin2d() +
    scale_fill_gradientn(colours=colvec) +
    facet_wrap(. ~ asm, ncol=1) +
    ggtitle('Telomereic Gene Hits') +
    scale_x_continuous(name='Number of Hits', limits=c(0,10), breaks=seq(0,10,1)) +
    scale_y_continuous(name='Adjusted Copy Number', limits=c(0,10), breaks=seq(0,10,1)) +
    theme_bw()
print(asmheatmap)
dev.off()
copynumcorrpdf=file.path(dbxdir, 'paperfigs', 'raw', 'copynum_corr.pdf')
pdf(copynumcorrpdf, height=5, width=9)
corrheat=ggplot(allhitscorr, aes(x=asm, y=ref)) +
    geom_bin2d() +
    scale_fill_gradientn(colours=colvec) +
    ggtitle('Adjusted Copy Number') +
    scale_x_continuous(name='JHU_Cniv_v1', limits=c(0,11), breaks=seq(0,11,1)) +
    scale_y_continuous(name='C. Nivariensis Reference', limits=c(0,11), breaks=seq(0,11,1)) +
    theme_bw()
print(corrheat)
dev.off()





repgenes=regiontable %>%
    filter(Is_Repetitive=='yes')
foundreps=intersect(repgenes$Systematic_Name, unique(nivhits$gene))

##filter gff for longer genes, more likely to be gpicwps
cnames=c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
annot=read_tsv(asmgff, col_names=cnames) %>%
    mutate(width=end-start) %>%
    filter(width>1000)
longpath=file.path(datadir, 'patho', 'long_genes.gff')
write_tsv(annot, longpath, col_names=FALSE)


##long_genes.fa made using gffread
aafile=file.path(datadir, 'patho', 'long_genes.fa')
aa=readAAStringSet(aafile,  format='fasta')
for (i in 1:34) {
    file=file.path(datadir, 'patho', 'predgpi_lists', paste0('long_genes', as.character(i), '.fa'))
    first=100*(i-1)+1
    last=min(100*i, length(aa))
    set=aa[first:last]
    writeXStringSet(set, file, format='fasta')
}





##get nucleotide fasta of relevant genes
gpifile=file.path(datadir, 'patho', 'predgpi.fasta')
gpiaa=readAAStringSet(gpifile, format='fasta')
cnames=c('transcript', 'fp', 'omega')
gpilist=as_tibble(matrix(unlist(strsplit(names(gpiaa), split=' | ', fixed=TRUE)), ncol=3, byrow=TRUE))
colnames(gpilist)=cnames

cds=annot %>%
    filter(feature=='CDS')

regions_from_name <- function(name) {
    ##add genomic loci to a gpilist row, duplicate rows if necessary
    cdsinfo=cds[grepl(name$transcript, cds$attribute, fixed=TRUE),] %>%
        mutate(fp=name$fp) %>%
        mutate(omega=name$omega)
    return(cdsinfo)
}
    
annotgpi=gpilist %>%
    group_by(transcript) %>%
    do(regions_from_name(.))
asm=readDNAStringSet(asmfile)

gpiregions=GRanges(seqnames=annotgpi$chr,
                   ranges=IRanges(start=annotgpi$start,
                                  end=annotgpi$end,
                                  names=annotgpi$transcript),
                   strand=annotgpi$strand)
gpiseqs=getSeq(asm, gpiregions)
gpiseqfile=file.path(datadir, 'patho', 'gpigenes.fasta')
writeXStringSet(gpiseqs, gpiseqfile)




##read blast of gpi and compare
nivpredtsv=file.path(datadir, 'patho', 'nivar.final.predgpi_hits.tsv')
refpredtsv=file.path(datadir, 'patho', 'candida_nivariensis.predgpi_hits.tsv')
cnames=c('gene', 'chr', 'ident', 'alignlen', 'mismatches', 'gap', 'qstart', 'qend', 'sstart', 'send', 'eval', 'store')
nivpred=read_tsv(nivpredtsv, comment='#', col_names=cnames) %>%
    mutate(genome='asm') %>%
    rowwise() %>%
    mutate(start=min(sstart, send)) %>%
    mutate(end=max(sstart, send))
refpred=read_tsv(refpredtsv, comment='#', col_names=cnames) %>%
    mutate(genome='ref')
allpred=rbind(nivpred, refpred) %>%
    group_by(gene) %>%
    summarise(maxref=max(alignlen[genome=='ref']),
              maxasm=max(alignlen[genome=='asm']),
              numref=sum(genome=='ref'),
              numasm=sum(genome=='asm'))



##get list of pred genes that have repeats
trffile=file.path(datadir, 'trf', 'trf.out')
trfinfo=read_csv(trffile, col_names=c('trf'))
cnames=c('start', 'end', 'size', 'copies', 'consensus', 'matches', 'indels', 'score', 'A', 'C', 'G', 'T', 'entropy', 'seq1', 'seq2', 'seq3', 'seq4', 'name')
name=substr(trfinfo$trf[1], 2, str_length(trfinfo$trf[1]))
trf=as_tibble(str_split_fixed(trfinfo$trf[2], pattern=' ', n=17)) %>%
    mutate(name=name)
colnames(trf)=cnames
for (i in 3:dim(trfinfo)[1]) {
    if (substr(trfinfo$trf[i], 1, 1)=='@') {
        name=substr(trfinfo$trf[i], 2, str_length(trfinfo$trf[i]))
    } else {
        trfsep=as_tibble(str_split_fixed(trfinfo$trf[i], pattern=' ', n=17)) %>%
            mutate(name=name)
        colnames(trfsep)=cnames
        trf=rbind(trf, trfsep)
    }
}

repregs=GRanges(seqnames=trf$name,
                ranges=IRanges(start=as.numeric(trf$start), end=as.numeric(trf$end)))
predregs=GRanges(seqnames=nivpred$chr,
                 ranges=IRanges(start=nivpred$start, end=nivpred$end))
predrep=unique(nivpred$gene[unique(queryHits(findOverlaps(predregs, repregs)))])

repgpis=allpred %>%
    filter(gene %in% predrep)

gpipdf=file.path(dbxdir, 'paperfigs', 'raw', 'gpis.pdf')
pdf(gpipdf, h=6, w=18)
alen=ggplot(repgpis, aes(x=maxasm, y=maxref)) +
    geom_bin2d(binwidth=c(100,100)) +
    scale_fill_gradientn(colours=colvec) +
    ggtitle('Max GPI-CWP Hit Lengths') +
    scale_x_continuous(name='JHU_Cniv_v1', breaks=seq(0,12000,1000)) +
    scale_y_continuous(name='C. Nivariensis Reference', breaks=seq(0,4000,1000)) +
    theme_bw()
print(alen)
numalign=ggplot(repgpis, aes(x=numasm, y=numref)) +
    geom_bin2d(binwidth=c(100,100)) +
    scale_fill_gradientn(colours=colvec) +
    ggtitle('Max GPI-CWP Hit Lengths') +
    scale_x_continuous(name='JHU_Cniv_v1', breaks=seq(0,1300,100)) +
    scale_y_continuous(name='C. Nivariensis Reference', breaks=seq(0,100,100)) +
    theme_bw()
print(numalign)
dev.off()
    
