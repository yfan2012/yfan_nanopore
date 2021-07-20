#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr
dbxdir=~/Dropbox/timplab_data/mdr5
srcdir=~/Code/utils/qc

if [ $1 == check_flye ] ; then
    touch $dbxdir/qc/flye_asms.csv
    for i in MDRstool_16 MDRstool_19 ;
    do
	for len in 100m 10g 150m 500m ;
	do
	    python $srcdir/asm_assess.py \
		   -i $datadir/$i/metaflye/$len/${i}_${len}.assembly.fasta \
		   -p ${i}_${len} \
		   >> $dbxdir/qc/flye_asms.csv
	done
    done
fi


if [ $1 == check_run ] ; then

    for i in MDRstool_19 MDRstool_16 ;
    do
	Rscript $srcdir/run_summary.R -i $datadir/$i/called/sequencing_summary.txt -o $dbxdir/qc/${i}_run.pdf
    done
fi

    
if [ $1 == stool16 ] ; then

    touch $dbxdir/qc/stool16_asms.csv
    for i in native pcr ;
    do
	asm=$datadir/MDRstool_16/metaflye/$i/MDRstool_16_${i}.assembly.fasta
	python $srcdir/asm_assess.py \
	       -i $asm \
	       -p stool16_$i \
	       >> $dbxdir/qc/stool16_asms.csv
    done
    python $srcdir/asm_assess.py \
	   -i $datadir/MDRstool_16/metaspades/181127_hiC_stool_shotgun.contigs.fasta \
	   -p stool16_illumina \
	   >> $dbxdir/qc/stool16_asms.csv
fi
    
