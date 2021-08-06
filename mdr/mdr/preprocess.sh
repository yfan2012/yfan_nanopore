#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/mdr
prefix=200708_mdr_stool16native

ref=$datadir/ref/mdr_refs.fa
asm=$datadir/flye/$prefix/$prefix.assembly.fasta

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    cat $datadir/called/pass/*fastq.gz > $datadir/fastqs/$prefix.fq.gz
fi

if [ $1 == flye ] ; then
    mkdir -p $datadir/flye
    mkdir -p $datadir/flye/$prefix

    flye \
	--nano-raw $datadir/fastqs/$prefix.fq.gz \
	-o $datadir/flye/$prefix \
	-t 36 \
	-g 100m \
	--plasmids \
	--meta
    
    mv $datadir/flye/$prefix/assembly.fasta $datadir/flye/$prefix/$prefix.assembly.fasta
fi

if [ $1 == amr ] ; then
    mkdir -p $datadir/amr

    for i in ecoh card ncbi resfinder plasmidfinder vfdb ecoli_vf megares argannot ;
    do
	abricate \
	    --threads 36 \
	    --db $i \
	    $asm > $datadir/amr/$prefix.$i.tsv
    done
fi

spadesasm=/uru/Data/Nanopore/projects/mdr/MDRstool_16/metaspades/181127_hiC_stool_shotgun.contigs.fasta
if [ $1 == illumina_amr ] ; then
    prefix=`basename $spadesasm .contigs.fasta`
    for i in ecoh card ncbi resfinder plasmidfinder vfdb ecoli_vf megares argannot ;
    do
	abricate \
	    --threads 36 \
	    --db $i \
	    $spadesasm > $datadir/amr/$prefix.$i.tsv
    done
fi
