#!/bin/bash

datadir=/atium/Data/Nanopore/cpowgs
refdir=/mithril/Data/NGS/Reference/bcepacia
outdir=/atium/Data/Nanopore/cpowgs/refdbs
seqsdir=$datadir/References/gene_seqs

for i in $refdir/bcepacia*.fa ;
do
    ref=` echo $i | cut -d '/' -f 7 `
    strainraw=` echo $ref | cut -d '_' -f 2`
    strain=` echo $strainraw | cut -d '.' -f 1`

    makeblastdb -in $i -out $outdir/$strain.db -dbtype nucl

    blastn -query $seqsdir/all_genes.fa -db $outdir/$strain.db -outfmt 7 -out $datadir/genomes/var_profiles/$strain.genes_orig.tsv
    blastn -query $seqsdir/all_genes.fa -subject $i -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$strain.genes_orig.btop
    blastn -query $seqsdir/all_genes.fa -subject $i -dust no -parse_deflines -out $datadir/genomes/var_profiles/$strain.genes_orig.vis
    
done

    
