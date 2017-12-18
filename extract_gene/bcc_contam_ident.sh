#!/bin/bash

##make blast db from all penA genes in bcc
##blast hometown burkle against database
##win

seqsdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc/ident/contam_genome
dbdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refdb



##%identity between whole genome contam and others
contam=$seqsdir/bcontaminans.fa

for i in $seqsdir/*fa ;
do
    prefix=`echo ${i%.fa} | cut -d '/' -f 8 `
    
    blastn -query $i -db $dbdir/bcontaminans.db -outfmt 7 -out $outdir/$prefix.tsv
    blastn -query $i -subject $contam -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $outdir/$prefix.btop
    blastn -query $i -subject $contam -dust no -parse_deflines -out $outdir/$prefix.vis
done

