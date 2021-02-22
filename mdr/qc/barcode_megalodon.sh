#!/bin/bash

datadir=~/data/mdr/qc/megalodon
samps='neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm'

if [ $1 == index_megalodon_output ] ; then
    for i in $samps ;
    do
	time (python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
		-i $datadir/$i/$i/${i}_mod_basecalls.txt \
		-o $datadir/$i/$i/${i}_mod_basecalls.txt.idx ) \
	    &> $datadir/$i/$i/${i}_mod_index_time.txt
    done
fi


		
	     
