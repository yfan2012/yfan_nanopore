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

if [ $1 == countsnps ] ; then
    Rscript ~/Code/utils/count_snps.R -i $outdir/parsnp/parsnp.vcf -o $outdir/parsnp
fi

if [ $1 == canuparnsp ] ; then
    ##just want to give nina some basic numbers. including wtdbg2 tigs rn could be a whole *thing*
    mkdir - $outdir/parnsp_canu
    ref=$datadir/assemblies/parsnp/st5317_canu.pilon.29.fasta
    parsnp -r $ref -d $datadir/assemblies/parsnp -p 36 -o $outdir/parsnp_canu -c
    harvesttools -i $outdir/parsnp_canu/parsnp.ggr -V $outdir/parsnp_canu/parsnp.vcf
    Rscript ~/Code/utils/count_snps.R -i $outdir/parsnp_canu/parsnp.vcf -o $outdir/parsnp_canu
fi


