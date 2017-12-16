#!/bin/bash

datadir=/atium/Data/Nanopore/cpowgs/170816_BUCC
illdir=$datadir/bcc_align_ill
nanodir=$datadir/bcc_align

for i in $illdir/*.sorted.bam ;
do
    (prefix=${i%.sorted.bam}
    samtools flagstat $i > $prefix.ill.flagstat.txt) &
done

for i in $nanodir/*.sorted.bam ;
do
    (prefix=${i%.sorted.bam}
    samtools flagstat $i > $prefix.nano.flagstat.txt ) &
done
