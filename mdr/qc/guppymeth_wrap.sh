#!/bin/bash

samps='neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm'

if [ $1 == test ] ; then
    time ( bash megalodon.sh nebdcm rerio_allmod ) &> test_time.txt
    ##bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm CCWGG 1 megalodon_vanilla
fi

if [ $1 == rerio ] ; then
    for i in $samps ;
    do
	echo $i
	time ( bash megalodon.sh $i rerio_allmod ) &> ~/data/mdr/qc/megalodon/$i.time.txt
    done
fi

	
