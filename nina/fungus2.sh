#!/bin/bash

srcdir=~/Code/utils/marcc
datadir=/scratch/groups/mschatz1/cpowgs/fungus/


if [ $1 == untar ] ; then
    ##NB manually moved stuff to mike's storage allocation. These paths may be wrong
    mkdir $datadir/raw
    mkdir $datadir/batch_logs
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_nina $srcdir/untar.scr $datadir/180827_nina_fungus2.tar.gz $datadir
fi

if [ $1 == untar_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2
    mkdir -p $datadir/181108_nina_v2/raw
    mkdir -p $datadir/181108_nina_v2/batch_logs
    ##sbatch --output=$datadir/181108_nina_v2/batch_logs/untar.out --job-name=ut_nina $srcdir/untar.scr $datadir/181108_nina_v2.tar.gz $datadir/181108_nina_v2
    ##shitty lustre file system caused untar to time out
    bash $srcdir/untar.scr $datadir/181108_nina_v2.tar.gz $datadir/181108_nina_v2
fi
    
if [ $1 == call ] ; then
    mkdir $datadir/called
    mkdir $datadir/call_logs
    mkdir $datadir/call_done
    sbatch --array=0-141 --job-name=nina_call --output=$datadir/call_logs/nina_call.%A_%a.out $srcdir/bc_call.scr $datadir
fi

if [ $1 == call_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/called
    mkdir -p $datadir/181108_nina_v2/call_logs
    mkdir -p $datadir/181108_nina_v2/call_done
    sbatch --array=0-1901 --job-name=nina_call --output=$datadir/181108_nina_v2/call_logs/nina_call.%A_%a.out $srcdir/bc_call_LSK109.scr $datadir/181108_nina_v2
fi

if [ $1 == fastq ] ; then
    mkdir $datadir/fastqs
    cat $datadir/called/*/workspace/pass/barcode04/*.fastq > $datadir/fastqs/st90853.fastq
    cat $datadir/called/*/workspace/pass/barcode05/*.fastq > $datadir/fastqs/st31.fastq
fi

if [ $1 == fastq_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/fastqs/
    cat $datadir/181108_nina_v2/called/*/workspace/pass/barcode12/*.fastq > $datadir/181108_nina_v2/fastqs/st90853.fastq
    cat $datadir/181108_nina_v2/called/*/workspace/pass/barcode11/*.fastq > $datadir/181108_nina_v2/fastqs/st31.fastq
fi



if [ $1 == assemble_90853_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/st90853_assembly
    canu \
	-p st90853 -d $datadir/181108_nina_v2/st90853_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=39m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/181108_nina_v2/fastqs/st90853.fastq
fi


if [ $1 == assemble_90853 ] ; then
    mkdir -p $datadir/st90853_assembly
    ~/software/canu-1.7/Linux-amd64/bin/canu \
	-p st90853 -d $datadir/assembly \
	-gridOptions="--mem=8g --time=22:00:00 --account=mschatz1" \
	-gnuplotTested=True \
	-minReadLength=500 \
	genomeSize=39m \
	-nanopore-raw $datadir/fastqs/st90853.fastq
fi

if [ $1 == assemble_31_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/st31_assembly
    canu \
	-p st31 -d $datadir/181108_nina_v2/st31_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=39m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/181108_nina_v2/fastqs/st31.fastq
fi

if [ $1 == assemble_31 ] ; then
    mkdir -p $datadir/st31_assembly
    ~/software/canu-1.7/Linux-amd64/bin/canu \
	-p st31 -d $datadir/st31_assembly \
	-gridOptions="--mem=8g --time=22:00:00 --account=mschatz1" \
	-gnuplotTested=True \
	-minReadLength=500 \
	genomeSize=39m \
	-nanopore-raw $datadir/fastqs/st31.fastq
fi
