#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/fungus/1704XX_Lprolif

if [ $1 == np_long ] ; then
    python ~/Code/utils/fastq_long.py -i $datadir/fastqs/Lprolif_npreads.fastq -o $datadir/fastqs/Lprolif_npreads_over3kb.fastq -l 3000
fi

if [ $1 == canu ] ; then
    mkdir -p $datadir/canu_assembly

    canu \
	-p st5317 -d $datadir/canu_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=39m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/fastqs/Lprolif_npreads_over3kb.fastq
fi

if [ $1 == trim_illumina ] ; then
    java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 \
	 $datadir/fastqs/170417-scedo-hyphae_S1_L001_R1_001.fastq.gz \
	 $datadir/fastqs/170417-scedo-hyphae_S1_L001_R2_001.fastq.gz \
	 $datadir/fastqs/scedo_forward_paired.fq.gz \
	 $datadir/fastqs/scedo_forward_unpaired.fq.gz \
	 $datadir/fastqs/scedo_reverse_paired.fq.gz \
	 $datadir/fastqs/scedo_reverse_unpaired.fq.gz \
	 ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
fi

if [ $1 == wtdbg2 ] ; then
    mkdir -p $datadir/wtdbg2_assembly
    wtdbg2 -t 36 -i $datadir/fastqs/Lprolif_npreads_over3kb.fastq -fo $datadir/wtdbg2_assembly/Lprolif_wtdbg2
    wtpoa-cns -t 36 -i $datadir/wtdbg2_assembly/Lprolif_wtdbg2.ctg.lay.gz -fo $datadir/wtdbg2_assembly/Lprolif.wtdbg2.contigs.fasta
fi

if [ $1 == pilon ] ; then
    ##for i in canu wtdbg2 ;
    for i in wtdbg2 canu ;
    do
	mkdir -p $datadir/${i}_pilon
	##cp $datadir/fastqs/*paired* $datadir/${i}_pilon/
	sbatch pilon.scr $datadir/${i}_pilon $datadir/${i}_assembly/*.contigs.fasta st5317 $i
    done
fi

	

	
