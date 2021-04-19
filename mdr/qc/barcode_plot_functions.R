library(tidyverse)
library(RColorBrewer)
library(cowplot)

plot_bc_dists <- function(samp){
    ##takes samp info and plots dists
    bcfile=file.path(datadir, paste0(samp$samp, '_barcodes.txt'))
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


