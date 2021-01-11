library(ggplot2)
library(tidyverse)
library(gridExtra)

##quick pie chart for ubiome talk
##csv taken from pavian output
infected=gather(read_csv('~/Dropbox/yfan/griffin_sindbis/sindbis_infected_pavian.csv'))
pie=infected[c(6,8,9),]
pie=add_row(pie, key='Other', value=as.character(100-sum(as.numeric(pie$value))))

pie=pie %>%
    mutate(prop=as.numeric(value))
 
chart=ggplot(pie, aes(x='', y=prop, fill=key)) +
    geom_bar(width=1, stat='identity') +
    coord_polar('y', start=0) +
    ggtitle('Read classification' ) +
    theme_bw()

pdf('~/Dropbox/yfan/griffin_sindbis/infected_pie.pdf', width=8, height=8)
print(chart)
dev.off()

pdf('~/Dropbox/yfan/griffin_sindbis/infected_table.pdf', width=8, height=8)
print(grid.table(pie[,c(1,3)], rows=NULL))
dev.off()
