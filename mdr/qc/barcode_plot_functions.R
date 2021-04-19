library(tidyverse)
library(RColorBrewer)
library(cowplot)

plot_bc_dists <- function(bcfile){
    ##takes samp info and plots dists
    bc=read_tsv(bcfile, col_names=bc_cols) %>%
        mutate(samp=samp$motif) %>%
        gather('motif', 'count', -readname, -samp)
    
    plot=ggplot(bc, aes(x=count, colour=motif, fill=motif, alpha=.2)) +
        geom_density() +
        ggtitle(paste0(samp$samp, '_', samp$motif)) +
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
    
    
