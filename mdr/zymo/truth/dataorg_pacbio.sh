#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/zymo/truth

if [ $1 == getinfo ] ; then

    awk -F , '{print $1, $12}' $datadir/pacbio/raw/sra_metadata.csv > $datadir/pacbio/raw/pacbio_names.txt

fi

if [ $1 == rename ] ; then
    while read p ; do
	acc=`echo $p | cut -d ' ' -f 1`
	name=`echo $p | cut -d ' ' -f 2`
	if [ $acc != Run ] ; then
	    mv $datadir/pacbio/raw/$acc $datadir/pacbio/raw/$name
	fi
    done < $datadir/pacbio/raw/pacbio_names.txt

fi

if [ $1 == get_pb_labels ] ; then
    Rscript ~/Code/yfan_nanopore/mdr/zymo/truth/pacbio_labels.R
fi


if [ $1 == cp_modinfo ] ; then
    mkdir -p $datadir/pacbio/modinfo
    
    while read p ; do
	number=`echo $p | cut -d ' ' -f 1 `
	label=`echo $p | cut -d ' ' -f 2 `

	cp $datadir/pacbio/smrtanalysis/016/$number/data/modifications.csv.gz $datadir/pacbio/modinfo/$label.modifications.csv.gz
	cp $datadir/pacbio/smrtanalysis/016/$number/data/motifs.gff.gz $datadir/pacbio/modinfo/$label.motifs.gff.gz
	cp $datadir/pacbio/smrtanalysis/016/$number/data/motif_summary.csv $datadir/pacbio/modinfo/$label.motif_summary.csv
    done < ~/Code/yfan_nanopore/mdr/zymo/truth/pblist_withlabels.txt
fi

	
    
