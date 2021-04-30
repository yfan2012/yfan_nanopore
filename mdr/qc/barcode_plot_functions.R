library(tidyverse)
library(RColorBrewer)
library(cowplot)

plot_bc_dists <- function(bcfile, title){
    ##takes samp info and plots dists
    bc_cols=c('readname', 'GATC', 'GANTC', 'CCWGG', 'GCNGC')
    bc=read_tsv(bcfile, col_names=bc_cols) %>%
        mutate(samp=samp$motif) %>%
        gather('motif', 'count', -readname, -samp)
    
    plot=ggplot(bc, aes(x=count, colour=motif, fill=motif, alpha=.2)) +
        geom_density() +
        ggtitle(title) +
        xlab('Barcode Score') +
        scale_x_log10() +
        scale_fill_brewer(palette = 'Set2') +
        scale_colour_brewer(palette = 'Set2') +
        theme_bw()
    return(plot)
}

plot_bc_dists_fromtibble <- function(bc, title){
    bc=bc %>%
        mutate(samp=samp$motif) %>%
        gather('motif', 'count', -readname, -samp)
    
    plot=ggplot(bc, aes(x=count, colour=motif, fill=motif, alpha=.2)) +
        geom_density() +
        ggtitle(title) +
        xlab('Barcode Score') +
        scale_x_log10() +
        scale_fill_brewer(palette = 'Set2') +
        scale_colour_brewer(palette = 'Set2') +
        theme_bw()
    return(plot)
}


plot_barnyard <- function(bcfile1, bcfile2,  motif1, motif2) {
    ##takes two samps and plots barnyard
    bc_cols=c('readname', 'GATC', 'GANTC', 'CCWGG', 'GCNGC')
    bc1=read_tsv(bcfile1, col_names=bc_cols) %>%
        mutate(samp=motif1)
    bc2=read_tsv(bcfile2, col_names=bc_cols) %>%
        mutate(samp=motif2)
    bcinfo=bind_rows(bc1, bc2)
    
    
    plot=ggplot(bcinfo, aes(x=get(motif1), y=get(motif2), colour=samp, alpha=.1)) +
        geom_point() +
        xlab(motif1) +
        ylab(motif2) +
        scale_colour_brewer(palette='Set2') +
        theme_bw()
}
    
call_part_avg <- function(bcinfo, thresh) {
    ##thresh is expressed as +/- proportion of average
    ##output [readname, chrname, bin_barcode]
    scores=bcinfo[,-(1:2)]
    minscores=rowMeans(scores)*thresh
    test=scores>minscores
    barcode=apply(test, 1, function(x) paste0(as.character(as.numeric(x)), collapse=''))

    barcodedreads=bcinfo %>%
        select(readname, chrname) %>%
        mutate(barcode=barcode)

    return(barcodedreads)
}


plot_pops <- function(pops, name) {
    plot=ggplot(pops, aes(x=barcode, y=pops, colour=samp, fill=samp, alpha=.2)) +
        geom_bar(stat='identity') +
        ggtitle(name) +
        scale_color_brewer(palette='Set2') +
        scale_fill_brewer(palette='Set2') +
        theme_bw()
    return(plot)
}
