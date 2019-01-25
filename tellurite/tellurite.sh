#!/bin/bash

illrawdir=/work-zfs/mschatz1/cpowgs/illumina_tellurite

if [ $1 == dl_illumina ] ; then
    fastq-dump -O $illrawdir/ -A SRR7128305
    mv $illrawdir/SRR7128305.fastq $illrawdir/kp590.fastq

    fastq-dump -O $illrawdir/ -A SRR7128294
    mv $illrawdir/SRR7128294.fastq $illrawdir/kp697.fastq

    fastq-dump -O $illrawdir/ -A SRR7128274
    mv $illrawdir/SRR7128274.fastq $illrawdir/kp723.fastq
    
    fastq-dump -O $illrawdir/ -A SRR7128297
    mv $illrawdir/SRR7128297.fastq $illrawdir/kp725.fastq

    fastq-dump -O $illrawdir/ -A SRR7128299
    mv $illrawdir/SRR7128299.fastq $illrawdir/kp860.fastq

    fastq-dump -O $illrawdir/ -A SRR7128280
    mv $illrawdir/SRR7128280.fastq $illrawdir/kp1165.fastq

    fastq-dump -O $illrawdir/ -A SRR7128341
    mv $illrawdir/SRR7128341.fastq $illrawdir/kp1469.fastq

    fastq-dump -O $illrawdir/ -A SRR7128348
    mv $illrawdir/SRR7128348.fastq $illrawdir/kp1497.fastq

    fastq-dump -O $illrawdir/ -A SRR7128339
    mv $illrawdir/SRR7128339.fastq $illrawdir/kp1559.fastq

    fastq-dump -O $illrawdir/ -A SRR7128284
    mv $illrawdir/SRR7128284.fastq $illrawdir/kp1606.fastq

    fastq-dump -O $illrawdir/ -A SRR7128363
    mv $illrawdir/SRR7128363.fastq $illrawdir/kp1637.fastq

    fastq-dump -O $illrawdir/ -A SRR7128365
    mv $illrawdir/SRR7128365.fastq $illrawdir/kp1693.fastq

    fastq-dump -O $illrawdir/ -A SRR7128285
    mv $illrawdir/SRR7128285.fastq $illrawdir/kp1946.fastq

    fastq-dump -O $illrawdir/ -A SRR7128328
    mv $illrawdir/SRR7128328.fastq $illrawdir/kp1958.fastq

    fastq-dump -O $illrawdir/ -A SRR7128352
    mv $illrawdir/SRR7128352.fastq $illrawdir/kp2007.fastq

    fastq-dump -O $illrawdir/ -A SRR7128290
    mv $illrawdir/SRR7128290.fastq $illrawdir/kp2026.fastq

    fastq-dump -O $illrawdir/ -A SRR712831
    mv $illrawdir/SRR7128313.fastq $illrawdir/kp2156.fastq
fi

if [ $1 == troubleshoot_dl ] ; then
    fastq-dump -O $illrawdir/ -A SRR7128313
    mv $illrawdir/SRR7128313.fastq $illrawdir/KLPN_2156.fastq
fi

datadir=/scratch/groups/mschatz1/cpowgs/tellurite
if [ $1 == circlator ] ; then
    for i in $datadir/KLPN* ;
    do
	
