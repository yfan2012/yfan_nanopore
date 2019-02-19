#!/bin/bash

srcdir=~/Code/utils/marcc
datadir=/scratch/groups/mschatz1/cpowgs/fungus/181107_birdpair


if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    mkdir -p $datadir/batch_logs
    ##sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_bird $srcdir/untar.scr $datadir/181107_birdpair.tar.gz $datadir
    bash $srcdir/untar.scr $datadir/181107_birdpair.tar.gz $datadir
    
fi

if [ $1 == call ] ; then
    mkdir $datadir/called
    mkdir $datadir/call_logs
    mkdir $datadir/call_done
    sbatch --array=0-1356 --job-name=birdpair --output=$datadir/call_logs/bird_call.%A_%a.out $srcdir/bc_call_LSK109.scr $datadir
fi

if [ $1 == fastq ] ; then
    mkdir -p $datadir/fastqs
    cat $datadir/called/*/workspace/pass/barcode09/*.fastq > $datadir/fastqs/patient.fastq
    cat $datadir/called/*/workspace/pass/barcode10/*.fastq > $datadir/fastqs/bird.fastq
fi

if [ $1 == long_fq ] ; then
    python ~/Code/utils/fastq_long.py -i $datadir/fastqs/patient.fastq -o $datadir/fastqs/patient_over5kb.fastq -l 5000
    python ~/Code/utils/fastq_long.py -i $datadir/fastqs/bird.fastq -o $datadir/fastqs/bird_over5kb.fastq -l 5000
fi

if [ $1 == wtdbg2 ] ; then
    mkdir -p $datadir/patient_wtdbg2
    mkdir -p $datadir/bird_wtdbg2
    wtdbg2 -t 32 -i $datadir/fastqs/patient_over5kb.fastq -fo $datadir/patient_wtdbg2/patient_wtdbg2
    wtdbg2 -t 32 -i $datadir/fastqs/bird_over5kb.fastq -fo $datadir/bird_wtdbg2/bird_wtdbg2
    wtpoa-cns -t 32 -i $datadir/patient_wtdbg2/patient_wtdbg2.ctg.lay.gz -fo $datadir/patient_wtdbg2/patient_wtdbg2.contigs.fasta
    wtpoa-cns -t 32 -i $datadir/bird_wtdbg2/bird_wtdbg2.ctg.lay.gz -fo $datadir/bird_wtdbg2/bird_wtdbg2.contigs.fasta
fi

if [ $1 == pilon ] ; then
    ##sed -i -e 's/ /_/g' $datadir/bird_wtdbg2/bird_wtdbg2.contigs.fasta
    ##sed -i -e 's/ /_/g' $datadir/patient_wtdbg2/patient_wtdbg2.contigs.fasta
    
    mkdir -p $datadir/pilon_patient
    mkdir -p $datadir/pilon_bird

    cp /work-zfs/mschatz1/cpowgs/fungus/illumina/patient* $datadir/pilon_patient/
    cp /work-zfs/mschatz1/cpowgs/fungus/illumina/bird* $datadir/pilon_bird/ 
    
    sbatch --output=$datadir/batch_logs/bird_pilon.out --job-name=bird ./pilon.scr $datadir/pilon_bird $datadir/bird_wtdbg2/bird_wtdbg2.contigs.fasta bird wtdbg2
    sbatch --output=$datadir/batch_logs/patient_pilon.out --job-name=patient ./pilon.scr $datadir/pilon_patient $datadir/patient_wtdbg2/patient_wtdbg2.contigs.fasta patient wtdbg2
fi
