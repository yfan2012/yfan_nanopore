#!/bin/bash

##Nevermind, this is way too slow. Just do mummer. Duh. 

datadir=/atium/Data/Nanopore/cpowgs/170816_BUCC
refdir=$datadir/bcc_refs
clustdir=$datadir/bcc_clust
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc/ident

mkdir -p $clustdir

for i in $refdir/*fa ;
do
    prefix=` echo ${i%.fa} | cut -d '/' -f 8 `

    cp $i $clustdir/$prefix.fa

    sed -i -e "s/>/>$prefix/g" "$clustdir/$prefix.fa"
done


cat $clustdir/*.fa > $clustdir/all_genomes.fa

sed -i -e "s/A>/A\n>/g" "$clustdir/all_genomes.fa"
sed -i -e "s/C>/C\n>/g" "$clustdir/all_genomes.fa"
sed -i -e "s/G>/G\n>/g" "$clustdir/all_genomes.fa"
sed -i -e "s/T>/T\n>/g" "$clustdir/all_genomes.fa"

clustalo -i $clustdir/all_genomes.fa --full --percent-id --distmat-out=$outdir/clust.mat -o $outdir/clust_out.fa --seqtype=DNA --threads=12 --force
