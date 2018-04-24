#!/bin/bash

refdir=/mithril/Data/NGS/Reference/bcc

for i in $refdir/*fa;
do
    fa=`echo $i | cut -d '/' -f 7`
    prefix=${fa%.fa}
    species=${prefix#b}
    wget -O /mithril/Data/NGS/Reference/bcc/annotations/$species.gbk ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Burkholderia_$species/assembly_summary.txt
done
