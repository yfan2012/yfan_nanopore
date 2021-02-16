#!/bin/bash

##getting data organized for qc of megalodon/guppy meth calling.
##mostly stolen from yfan_meth/rerio
##subset stuff is already set up pretty nicely in mithril

rawdir=/mithril/Data/Nanopore/projects/methbin
datadir=~/data/mdr

if [ $1 == cp_multiraw_sub ] ; then
    mkdir -p $datadir/qc
    cp -r $rawdir/multiraw_sub $datadir/qc
fi

    
    
