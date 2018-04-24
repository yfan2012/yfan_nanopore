#!/bin/bash

datadir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc
refdir=/atium/Data/Nanopore/cpowgs/References/ref_gene_seqs


for i in gyrA parC;
do
    transeq $datadir/$i.fasta $datadir/$i.forward.faa -frame=F
    transeq $datadir/$i.fasta $datadir/$i.reverse.faa -frame=R
    transeq $refdir/bcepacia_$i.fasta $datadir/$i.ref.forward.faa -frame=F
    transeq $refdir/bcepacia_$i.fasta $datadir/$i.ref.reverse.faa -frame=R
done


##Committed the sin of manual file manipulation.
blastp -query $datadir/ref_gyrA_for.faa -subject $datadir/tig2_gyrA_for.faa > $datadir/gyrA.vis
blastp -query $datadir/ref_parC_for.faa -subject $datadir/tig2_parC_for.faa > $datadir/parC_for.vis
blastp -query $datadir/ref_parC2_rev.faa -subject $datadir/tig2_parC_rev.faa > $datadir/parC_rev.vis
