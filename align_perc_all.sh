#!/bin/bash

outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc
datadir=/atium/Data/Nanopore/cpowgs/170816_BUCC

rm $outdir/align_rates.csv
touch $outdir/align_rates.csv
for i in $datadir/bcc_align/*flagstat.txt ;
do
    echo $i
    python ~/Code/carbapenem_r21/bucc/align_perc.py -i $i >> $outdir/align_rates.csv
done

rm $outdir/align_rates_ill.csv
touch $outdir/align_rates_ill.csv
for i in $datadir/bcc_align_ill/*flagstat.txt ;
do
    echo $i
    python ~/Code/carbapenem_r21/bucc/align_perc.py -i $i >> $outdir/align_rates_ill.csv
done

	 
