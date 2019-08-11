#!/bin/bash

datadir=/kyber/Data/Nanopore/projects/nina
asmdir=$datadir/assemblies/pilon
outdir=~/Dropbox/yfan/nina_fungus

if [ $1 == rename ] ; then
    sed -i -e 's/_pilon//g' $asmdir/*fasta
    for i in $asmdir/*fasta ;
    do
	stid=`basename $i .fasta`
	sed -i -e "s/>/>$stid/g" $i
    done
fi
    
if [ $1 == parsnp ] ; then
    ref=$asmdir/st5317_canu.pilon.29.fasta
    mkdir -p $outdir/parsnp
    parsnp -r $ref -d $asmdir -p 36 -o $outdir/parsnp -c
    harvesttools -i $outdir/parsnp/parsnp.ggr -V $outdir/parsnp/parsnp.vcf
fi

