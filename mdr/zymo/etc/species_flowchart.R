library(DiagrammeR)
library(DiagrammeRsvg)  # for conversion to svg
library(rsvg)
library(tidyverse)
library(RColorBrewer)

##draw species tree

dbxdir='~/gdrive/mdr/zymo'

mycolors=brewer.pal(8,'Set2')

speciespdf=file.path(dbxdir, 'species_chart.pdf')

chart=grViz(diagram = "digraph dot {
  graph[layout = dot, fontsize = 15]
  node [fontname = arial, shape = box]
  tab1 [label = '@@1']
  tab2 [label = '@@2']
  tab3 [label = '@@3']
  tab4 [label = '@@4']
  tab5 [label = '@@5']
  tab6 [style=filled, label = '@@6', fillcolor='@@22']
  tab7 [style=filled, label = '@@7', fillcolor='@@28']
  tab8 [label = '@@8']
  tab9 [style=filled, label = '@@9', fillcolor='@@25']
  tab10 [label = '@@10']
  tab11 [label = '@@11']
  tab12 [style=filled, label = '@@12', fillcolor='@@23']
  tab13 [label = '@@13']
  tab14 [label = '@@14']
  tab15 [label = '@@15']
  tab16 [label = '@@16']
  tab17 [style=filled, label = '@@17', fillcolor='@@24']
  tab18 [style=filled, label = '@@18', fillcolor='@@27']
  tab19 [label = '@@19']
  tab20 [label = '@@20']
  tab21 [style=filled, label = '@@21', fillcolor='@@26']

  tab1 -> tab2 -> tab3 -> tab4 -> tab5 -> tab6;
  tab5 -> tab7;
  tab4 -> tab8 -> tab9;
  tab3 -> tab10 -> tab11 -> tab12;
  tab1 -> tab13 -> tab14 -> tab15 -> tab16 -> tab17;
  tab16 -> tab18;
  tab14 -> tab19 -> tab20 -> tab21
}
  
  [1]: 'Bacteria'
  [2]: 'Firmicutes'    
  [3]: 'Bacilli'
  [4]: 'Bacillales'
  [5]: 'Bacillaceae'
  [6]: 'Bacillus Subtilis'
  [7]: 'Staphylococcus aureus'
  [8]: 'Listeriaceae'
  [9]: 'Listeria monocytogenes'
  [10]: 'Lactobacillales'
  [11]: 'Enterococcaceae'
  [12]: 'Enterococcaceae faecalis'
  [13]: 'Proteobacteria'
  [14]: 'Gammaproteobacteria'
  [15]: 'Enterobacterales'
  [16]: 'Enterobacteriaceae'
  [17]: 'Escherichia coli'
  [18]: 'Salmonella enterica'
  [19]: 'Pseudomonadales'
  [20]: 'Pseudomonadaceae'
  [21]: 'Pseudomonas aeruginosa'
  [22]: mycolors[1]
  [23]: mycolors[2]
  [24]: mycolors[3]
  [25]: mycolors[4]
  [26]: mycolors[5]
  [27]: mycolors[6]
  [28]: mycolors[7]
  ")


chart %>%
    export_svg() %>%
    charToRaw %>% 
    rsvg_pdf(speciespdf)

