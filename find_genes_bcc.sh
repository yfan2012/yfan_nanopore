#!/bin/bash

datadir=/atium/Data/Nanopore/cpowgs
refdir=$datadir/170816_BUCC/bcc_refs
outdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refdb
seqsdir=$datadir/References/gene_seqs

for i in $refdir/*.fa ;
do
    ref=` echo $i | cut -d '/' -f 8 `
    strainraw=` echo $ref | cut -d '_' -f 2`
    strain=` echo $ref |cut -d '.' -f 1`

    makeblastdb -in $i -out $outdir/$strain.db -dbtype nucl

    blastn -query $seqsdir/all_genes.fa -db $outdir/$strain.db -outfmt 7 -out $datadir/170816_BUCC/bcc_vars/$strain.genes_orig.tsv
    blastn -query $seqsdir/all_genes.fa -subject $i -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/170816_BUCC/bcc_vars/$strain.genes_orig.btop
    blastn -query $seqsdir/all_genes.fa -subject $i -dust no -parse_deflines -out $datadir/170816_BUCC/bcc_vars/$strain.genes_orig.vis
    
done

    
