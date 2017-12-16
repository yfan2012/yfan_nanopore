#!/bin/bash

annotdir=/mithril/Data/NGS/Reference/bcc/annotations

##From prokka website https://github.com/tseemann/prokka#databases
prokka-genbank_to_fasta_db --format=gff $annotdir/*.gbk > $annotdir/Burkholderia.faa

cdhit -i $annotdir/Burkholderia.faa -o $annotdir/Burkholderia -T 0 -M 0 -g 1 -s 0.8 -c 0.9
makeblastdb -dbtype prot -in $annotdir/Burkholderia
mv $annotdir/Burkholderia.p* ~/software/prokka/db/genus/

prokout=~/Dropbox/Lab/carbapenem_r21/annotations/bucc
prokka --outdir $prokout --genus Burkholderia --usegenus --prefix bucc_prokka --cpus 12 --force /atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs/bcepacia_170816_BUCC_pilon.fa &> $prokout/$samp.prokka_log.txt
