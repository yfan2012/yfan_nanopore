#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/sindbis
srcdir=~/Code/utils/marcc

if [ $1 == untar ] ; then
   for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
   do
       mkdir -p $datadir/$i/raw
       mkdir -p $datadir/$i/batch_logs
       sbatch --output=$datadir/$i/batch_logs/untar.txt $srcdir/untar.scr $datadir/*$i*.tar.gz $datadir/$i
   done
fi

if [ $1 == call ] ; then
    ##for i in antibody mock infected ;
    for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	mkdir -p $datadir/$i/called
	mkdir -p $datadir/$i/call_done

	numdirs=`find $datadir/$i/raw/* -maxdepth 0 -type d | wc -l `
	dummy=1
	maxdir=`expr $numdirs - $dummy`
	echo $maxdir
	sbatch --array=0-$maxdir --output=$datadir/$i/batch_logs/call.txt --job-name=call_$i $srcdir/call_rna.scr $datadir/$i
    done
fi
