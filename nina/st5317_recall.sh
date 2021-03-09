#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans
##just see if re-calling st5217 data helps

i=st5317

if [ $1 == untar ] ; then
    mkdir -p ~/data/prolificans/$i
    tar xzf $datadir/raw/170412_scedo_hyphae_raw.tar.gz -C ~/data/prolificans/$i/
fi

if [ $1 == basecall ] ; then
    guppy_basecaller \
	--recursive \
	--compress_fastq \
	--kit SQK-LSK108 \
	--flowcell FLO-MIN106 \
	-i ~/data/prolificans/$i/raw \
	-s ~/data/prolificans/$i/called \
	-x 'cuda:0'
fi

if [ $1 == cp_illumina ] ; then
    mkdir -p $datadir/$i/reads
    cp -r $datadir/old_data_st5317/reads/illumina $datadir/$i/reads/
fi

if [ $1 == fqs ] ; then
    mkdir -p $datadir/$i/reads
    mkdir -p $datadir/$i/reads/ont

    cat ~/data/prolificans/$i/called/*fastq.gz \
	> $datadir/$i/reads/ont/$i.fastq.gz
fi

if [ $1 == longfq ] ; then
    python3 ~/Code/utils/fastq_long.py \
	    -i $datadir/$i/reads/ont/$i.fastq.gz \
	    -o $datadir/$i/reads/ont/${i}_long.fastq.gz \
	    -l 3000
fi

if [ $1 == canu ] ; then
    mkdir -p $datadir/$i/asm/canu
    canu \
	-p $i \
	-d $datadir/$i/asm/canu \
	genomeSize=39m \
	-nanopore-raw $datadir/$i/reads/ont/${i}_long.fastq.gz
fi


if [ $1 == flye ] ; then
    mkdir -p $datadir/$i/asm/flye
    flye \
	--nano-raw $datadir/$i/reads/ont/${i}_long.fastq.gz \
	-o $datadir/$i/asm/flye \
	-g 39m \
	-t 54
fi


if [ $1 == nameflye ] ; then
    cp $datadir/$i/asm/flye/assembly.fasta $datadir/$i/asm/flye/$i.ctgs.fasta
	
fi
	
