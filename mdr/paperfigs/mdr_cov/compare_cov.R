library(tidyverse)



####compare cov and max
maxbasedfile='~/gdrive/mdr/paperfigs/contig_level/plotallrands.csv'
maxbased=read_csv(maxbasedfile) %>%
    mutate(method='max')
covbasedfile='~/gdrive/mdr/paperfigs/contig_level_cov/plotallrands_cov.csv'
covbased=read_csv(covbasedfile) %>%
    mutate(method='cov')
randlabs=bind_rows(maxbased, covbased) %>%
    filter(label!='rando')


plotfile='~/gdrive/mdr/paperfigs/contig_level/clinical_contig_clusters_compare_maxcov.pdf'
pdf(plotfile, h=8, w=11)
plot=ggplot(randlabs, aes(x=seqtogether, y=seqpure, colour=method)) +
    scale_colour_brewer(palette='Set2') +
    geom_step() +
    ggtitle('max vs cov') +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw()
print(plot)
dev.off()
