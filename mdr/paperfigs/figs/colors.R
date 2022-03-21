library(tidyverse)
library(RColorBrewer)


projdir='/mithril/Data/Nanopore/projects/methbin'
prefix='200708_mdr_stool16native_perf'
datadir=file.path(projdir, 'paperfigs/contig_level')
dbxdir='~/gdrive/mdr/paperfigs/figs'


##chrom to bin info
chrombinsfile='/mithril/Data/Nanopore/projects/methbin/paperfigs/contig_level/tigs2bins.tsv'
chrombins=read_tsv(chrombinsfile)


####get tree from methylation distance
exvec=c('ATGCAT', 'GTCGAC', 'GANTC', 'GTWWAC', 'AAGCTT', 'CTCGAG', 'CTGCAG', 'CCGCGG')
methfreq=readRDS(file.path(datadir, 'clin_methfreq.rds'))
freqs=methfreq %>%
    rowwise() %>%
    filter(!motif %in% exvec) %>%
    spread(key=motif, value=freq)

nummotifs=length(table(methfreq$motif))
keepchroms=names(table(methfreq$chrom)[table(methfreq$chrom)==nummotifs])

methchroms=methfreq %>%
    filter(chrom %in% keepchroms)

chrominfo=methchroms %>%
    spread(key=motif, value=freq) %>%
    filter(chrom %in% chrombins$rname)

matchrominfo=as.matrix(chrominfo %>% select(-chrom))
rownames(matchrominfo)=chrominfo$chrom

labelinfo=tibble(label=chrominfo$chrom) %>%
    filter(label %in% chrombins$rname) %>%
    rowwise() %>%
    mutate(bins=chrombins$bin[chrombins$rname==label])


####ok i'm really in it now. *sobs*
##20 colors from https://sashamaps.net/docs/resources/20-colors/
##replace black and white with tyrian purple and #008080
colors=c('#e6194B',
         '#3cb44b',
         '#ffe119',
         '#4363d8',
         '#f58231',
         '#911eb4',
         '#42d4f4',
         '#f032e6',
         '#bfef45',
         '#fabed4',
         '#66B990',
         '#dcbeff',
         '#9A6324',
         '#aaffc3',
         '#808000',
         '#ffd8b1',
         '#000075',
         '#a9a9a9',
         '#630330',
         '#008080')

##blending algo from here https://github.com/javierbyte/colorblendjs
##interactive version: https://javier.xyz/cohesive-colors/
##hex to rbg from here https://stackoverflow.com/questions/43911071/r-hex-to-rgb-converter
colorblend <- function(c1,c2,intensity) {
    a=as.vector(col2rgb(c1))
    b=as.vector(col2rgb(c2))

    newcol=c()
    for (i in 1:3) {
        if (a[i]<128) {
            val=2*b[i]*a[i]/255
        }else{
            val=(255-2*(255-a[i])*(255-b[i])/255)
        }
        value=(val*intensity+(a[i]*(1-intensity)))
        newcol=c(newcol, as.integer(round(value)))
    }
    newhex=rgb2col(red=newcol[1], green=newcol[2], blue=newcol[3], alpha=FALSE)
    return(newhex)
}

tint='#FFBA70'
intensity=.3
newcolors=NULL
for (i in colors){
    newcolor=colorblend(i, tint, intensity)
    newcolors=c(newcolors, newcolor)
}

mycolors=tibble(bin=names(table(labelinfo$bins)),
                color=sample(colors))

write_csv(mycolors, 'colors.csv')








