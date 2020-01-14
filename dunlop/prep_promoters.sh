#!/bin/bash

refdir=/uru/Data/dunlop/reference

if [ $1 == preptxt ] ; then
    wget -O $refdir/regulondb_PromoterSet.txt http://regulondb.ccg.unam.mx/menu/download/datasets/files/PromoterSet.txt
    sed -i -e 's/TSS_/TSS/g' $refdir/regulondb_PromoterSet.txt
fi

if [ $1 == getfasta ] ; then
    ##get fasta of ecoli promoters from the regulondb promoterset file
    python annot_promoters.py -i $refdir/regulondb_PromoterSet.txt -o $refdir/regulondb_PromoterSet.fa
fi

if [ $1 == align_promoters ] ; then
    minimap2 $refdir/bw_BW25113.fa $refdir/regulondb_PromoterSet.fa > $refdir/bw_BW25113_promoters.paf
fi

if [ $1 == getgff ] ; then
    ##get gff of ecoli promoters from the alignment
    python annot_promoters.py -i $refdir/bw_BW25113_promoters.paf -o $refdir/bw_BW25113_promoters.gff
fi

if [ $1 == catgff ] ; then
    cat $refdir/bw_BW25113_promoters.gff > $refdir/bw_BW25113_withPromoters.gff3
    head --lines=-1 $refdir/bw_BW25113_promoters.gff >> $refdir/bw_BW25113_withPromoters.gff3
fi
