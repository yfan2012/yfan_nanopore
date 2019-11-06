#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar
dbxdir=~/Dropbox/yfan/nivar

if [ $1 == find_enrich ] ; then
    
    ref=$datadir/reference/candida_nivariensis.fa
    mkdir -p $dbxdir/motif_enrich
    
    for i in $datadir/mummer/*ref/*.snps ;
    do
	prefix=`basename $i .snps`
	echo $prefix
	python2 ~/Code/utils/motif_enrich.py -s $i -r $ref -m 6 -o $dbxdir/motif_enrich/$prefix.csv
    done

    for i in $datadir/mummer/*raw/*.snps ;
    do
	prefix=`basename $i .snps`
	echo $prefix
	corr=`echo $prefix | cut -d _ -f 3 `
	pore=`echo $prefix | cut -d _ -f 2`

	refcorr=$datadir/mummer/${pore}_${corr}_raw/$corr.15.fasta

	python2 ~/Code/utils/motif_enrich.py -s $i -r $refcorr -m 6 -o $dbxdir/motif_enrich/$prefix.csv
    done

fi

if [ $1 == venn_mummer ] ; then
    mkdir -p $datadir/mummer_venn

    for i in r9 r10
    do
	raw=$datadir/assemble/${i}_assembly/nivar_${i}.contigs.fasta

	rm ~/tmp/*

	for corr in pilon racon freebayes ;
	do

	    ##mummer against raw
	    mkdir -p $datadir/mummer_venn/${i}_raw_${corr}

	    cp $raw ~/tmp/
	    cp $datadir/$corr/${i}_$corr/*15*.fasta ~/tmp/$corr.15.raw.fasta
	    tr '[:lower:]' '[:upper:]' < ~/tmp/$corr.15.raw.fasta > ~/tmp/$corr.15.fasta
	    sed -i -e "s/>/>${i}_${corr}_15_/g" ~/tmp/$corr.15.fasta

	    nucmer -p ~/tmp/nivar_${i}_raw_${corr} ~/tmp/nivar_${i}.contigs.fasta ~/tmp/$corr.15.fasta 
	    mummerplot --filter --fat --png -p ~/tmp/nivar_${i}_raw_${corr} ~/tmp/nivar_${i}_raw_${corr}.delta
	    dnadiff -p ~/tmp/nivar_${i}_raw_${corr} ~/tmp/nivar_${i}.contigs.fasta ~/tmp/$corr.15.fasta 

	    cp ~/tmp/* $datadir/mummer_venn/${i}_raw_${corr}/

	    rm ~/tmp/*
	done
    done
fi

if [ $1 == error_venn ] ; then
    python2 ~/Code/utils/meth/error_venn.py \
	    -m1 $datadir/mummer/r9_raw_ref/nivar_r9_raw_ref.snps \
	    -m2 $datadir/mummer/r10_raw_ref/nivar_r10_raw_ref.snps \
	    -o $dbxdir/error_venn.csv
fi

if [ $1 == correction_venn ] ; then
    python2 ~/Code/yfan_nanopore/nivar/error_venn_3samps.py \
	    -m1 $datadir/mummer/r9_freebayes_ref/nivar_r9_freebayes_ref.snps \
	    -m2 $datadir/mummer/r9_pilon_ref/nivar_r9_pilon_ref.snps \
	    -m3 $datadir/mummer/r9_racon_ref/nivar_r9_racon_ref.snps \
	    -o $dbxdir/r9_correction_venn.csv
    
    python2 ~/Code/yfan_nanopore/nivar/error_venn_3samps.py \
	    -m1 $datadir/mummer/r10_freebayes_ref/nivar_r10_freebayes_ref.snps \
	    -m2 $datadir/mummer/r10_pilon_ref/nivar_r10_pilon_ref.snps \
	    -m3 $datadir/mummer/r10_racon_ref/nivar_r10_racon_ref.snps \
	    -o $dbxdir/r10_correction_venn.csv
fi
   

if [ $1 == neighbor ] ; then
    python2 ~/Code/yfan_nanopore/nivar/error_nearest.py \
	    -m1 $datadir/mummer/r9_raw_ref/nivar_r9_raw_ref.snps \
	    -m2 $datadir/mummer/r10_raw_ref/nivar_r10_raw_ref.snps \
	    -o $dbxdir/neighbor.csv
fi

if [ $1 == telos ] ; then
    r9=$datadir/freebayes/r9_freebayes/nivar_fb15_bwa.fasta
    r10=$datadir/freebayes/r10_freebayes/nivar_fb15_bwa.fasta
    grep -C 5 '>' $r9 > $dbxdir/telo_r9.txt
    grep -C 5 '>' $r10 > $dbxdir/telo_r10.txt
fi


if [ $1 == correction_types ] ; then
    touch $dbxdir/correction_types.csv
    rm $dbxdir/correction_types.csv
    touch $dbxdir/correction_types.csv

    echo pore,alg,deletions,insertions,snps >> $dbxdir/correction_types.csv
    
    for i in $datadir/mummer/*_raw/*snps ;
    do
	deletion=`awk '$3=="." {n++} END {printf("%d\n",n)}' $i`
	insertion=`awk '$2=="." {n++} END {printf("%d\n",n)}' $i`
	snps=`awk '$2!="." && $3!="." {n++} END {printf("%d\n",n)}' $i`

	pore=`basename $i .snps | cut -d _ -f2`
	alg=`basename $i .snps | cut -d _ -f3`
	
	echo $pore,$alg,$deletion,$insertion,$snps >> $dbxdir/correction_types.csv
    done
fi
