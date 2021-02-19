#!/bin/bash

datadir=~/data/mdr/qc/guppymeth
samps='neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm'

if [ $1 == call ] ; then
    mkdir -p $datadir
    for i in $samps ;
    do
	mkdir -p $datadir/$i/$i
	guppy_basecaller \
	    -i ~/data/mdr/qc/multiraw_sub/$i \
	    -s $datadir/$i/$i \
	    -d ~/software/rerio/basecall_models/ \
	    -c res_dna_r941_min_modbases-all-context_v001.cfg \
	    -x "cuda:0" \
	    --fast5_out
    done
fi
