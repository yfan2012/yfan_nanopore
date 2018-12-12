library(ggplot2)
library(tidyverse)

##read in cov files, basic cov plot
samps=c('samp15at24hrs','samp15at48hrs', 'samp20at24hrs', 'samp20at48hrs')
datadir='/scratch/groups/mschatz1/cpowgs/phage/'

allcov=NULL

for (samp in samps) {
    covfile=paste0(datadir, '/',samp, '/cov/', samp, '.primary.cov') 
    cov=read_tsv(covfile, col_names=c('chr', 'pos', 'cov')) %>%
        mutate(sample=samp)
    allcov=bind_rows(allcov, cov)
    
    xwidth=max(cov$pos)-min(cov$pos)
    yheight=max(cov$cov)
    
    pdf(paste0(datadir,'/', samp, '/cov/', samp, '.primary.cov.pdf'), width=20, height=5)
    print(ggplot(cov, aes(x=pos, y=cov)) +
        geom_line() +
        xlab('Sindbis Genome Position') +
        ylab('Depth') +
        ggtitle(paste0(samp, ' Coverage')) +
        theme_bw())
    dev.off()
}

pdf(paste0(datadir,'/allcov.primary.cov.pdf'), width=20, height=5)
print(ggplot(allcov, aes(x=pos, y=cov, colour=sample)) +
      geom_line() +
      xlab('Sindbis Genome Position') +
      ylab('Depth') +
      ggtitle(paste0('Coverage')) +
      theme_bw())
dev.off()
