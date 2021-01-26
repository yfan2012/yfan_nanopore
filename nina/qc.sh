#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans
dbxdir=~/Dropbox/timplab_data/prolificans

if [ $1 == mummer ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/mummer

	nucmer \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/genomes/$i.flye.fasta \
	    $datadir/$i/genomes/$i.canu.fasta
	
	mummerplot \
	    --filter --fat --postscript \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/mummer/$i.delta
	
	mummerplot \
	    --filter --fat --png \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/mummer/$i.delta

	dnadiff \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/genomes/$i.flye.fasta \
	    $datadir/$i/genomes/$i.canu.fasta
    done
fi

if [ $1 == asmstats ] ; then

    mkdir -p $dbxdir/qc
    touch $dbxdir/qc/asmstats.csv
    
    for i in st31 st90853 st5317 ;
    do
	python ~/Code/utils/qc/asm_assess.py \
	       -i $datadir/$i/genomes/$i.canu.fasta \
	       -p $i.canu >> $dbxdir/qc/asmstats.csv
	python ~/Code/utils/qc/asm_assess.py \
	       -i $datadir/$i/genomes/$i.flye.fasta \
	       -p $i.flye >> $dbxdir/qc/asmstats.csv
    done
fi


if [ $1 == runstats ] ; then

    touch $dbxdir/qc/runstats.csv
    touch $dbxdir/qc/runstats_long.csv
    
    for i in st31 st90853 st5317 ;
    do
	bash ~/Code/utils/qc/basic_run_assess.sh \
	     $datadir/$i/reads/ont/${i}.fastq.gz >> $dbxdir/qc/runstats.csv
	bash ~/Code/utils/qc/basic_run_assess.sh \
	     $datadir/$i/reads/ont/${i}_long.fastq.gz >> $dbxdir/qc/runstats_long.csv
    done
fi
