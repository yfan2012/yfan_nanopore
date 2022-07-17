#!/bin/bash

projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/paperfigs/nohuman
prefix=200708_mdr_stool16native_nohuman

if [ $1 == cp_nohumanraw ] ; then
    ##copy nohuman reads made during data cleaning to ssd drive
    cp -r $projdir/mdr/dataclean/$prefix ~/data/mdr/mdr/raw/
fi

if [ $1 == cp_nohumanfq ] ; then
    cp $projdir/mdr/dataclean/$prefix/$prefix.fq.gz $projdir/mdr/fastqs/$prefix.fq.gz
fi

    
