#!/bin/bash

datadir=~/work/180714_dunlop_3ecoli
srcdir=~/Code/utils/marcc

mkdir -p $datadir/batch_logs

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_dunlop $srcdir/untar.scr $datadir/180714_dunlop_3ecoli.tar.gz $datadir
fi

   
if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-693 --job-name=call_dunlop --output=$datadir/call_logs/dunlop_ecoli3.%A_%a.out $srcdir/bc_call.scr $datadir
fi


if [ $1 == fastq ] ; then
    mkdir -p $datadir/fastqs
    cat $datadir/called/*/workspace/pass/barcode01/*fastq > $datadir/fastqs/ecoli1.fastq
    cat $datadir/called/*/workspace/pass/barcode02/*fastq > $datadir/fastqs/ecoli2.fastq
    cat $datadir/called/*/workspace/pass/barcode03/*fastq > $datadir/fastqs/ecoli3.fastq
fi
    
if [ $1 == downsamp ] ; then
    head -800000 $datadir/fastqs/ecoli1.fastq > $datadir/fastqs/ecoli1_sub200k.fastq
    head -800000 $datadir/fastqs/ecoli2.fastq > $datadir/fastqs/ecoli2_sub200k.fastq
    head -800000 $datadir/fastqs/ecoli3.fastq > $datadir/fastqs/ecoli3_sub200k.fastq
fi

if [ $1 == assemble ] ; then
    mkdir $datadir/ecoli1_assembly
    bash assemble_bacteria.sh $datadir/fastqs/ecoli1.fastq $datadir/ecoli1_assembly
    mkdir $datadir/ecoli2_assembly
    bash assemble_bacteria.sh $datadir/fastqs/ecoli2.fastq $datadir/ecoli2_assembly
    mkdir $datadir/ecoli3_assembly
    bash assemble_bacteria.sh $datadir/fastqs/ecoli3.fastq $datadir/ecoli3_assembly
fi

if [ $1 == assemble17 ] ; then
    mkdir $datadir/ecoli1_assembly17
    bash assemble_bacteria.sh $datadir/fastqs/ecoli1_sub200k.fastq $datadir/ecoli1_assembly17
    mkdir $datadir/ecoli2_assembly17
    bash assemble_bacteria.sh $datadir/fastqs/ecoli2_sub200k.fastq $datadir/ecoli2_assembly17
    mkdir $datadir/ecoli3_assembly17
    bash assemble_bacteria.sh $datadir/fastqs/ecoli3_sub200k.fastq $datadir/ecoli3_assembly17
fi

    
if [ $1 == illumina_dilith ] ; then
    ##copy illumina data to dilithium
    mkdir -p /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli
    for i in `find ~/BaseSpace/Projects/180722_dunlop_3ecoli/ -name *fastq.gz` ;
    do
	cp $i /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/
    done
fi

if [ $1 == illumina_marcc ] ; then
    mkdir -p ~/work/180714_dunlop_3ecoli/illumina
    scp -r smaug:/dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/* ~/work/180714_dunlop_3ecoli/illumina/
fi
