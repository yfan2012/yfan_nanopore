#!/bin/bash

datadir=/kyber/Data/Nanopore/projects/tellurite

if [ $1 == mkdb ] ; then
    mkdir -p $datadir/reads_blast_db
    for i in $datadir/reads/*.fq ;
    do
	samp=`basename $i .fq`
	if [ ! -f $datadir/readsdb/$samp.fa ] ; then
	    echo $samp
	    seqtk seq -a $i > $datadir/reads_blast_db/$samp.fa
	    makeblastdb -in $datadir/reads_blast_db/$samp.fa -out $datadir/reads_blast_db/$samp.db -dbtype nucl
	fi
    done
fi

if [ $1 == ter ] ; then
    mkdir -p $datadir/reads_blast_hits
    cat $datadir/ref/*fa > $datadir/ref/klpn_genes.fasta
    for i in $datadir/reads/KLPN*.fq ;
    do
	samp=`basename $i .fq`
	if [ ! -f $datadir/blast_hits_reads/$samp.tsv ] ; then
	    echo $samp
	    blastn -query $datadir/ref/klpn_genes.fasta -db $datadir/reads_blast_db/$samp.db -outfmt 7 -out $datadir/reads_blast_hits/$samp.tsv
	fi
    done
fi
