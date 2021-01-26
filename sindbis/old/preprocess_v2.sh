#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis
srcdir=~/Code/utils/marcc

if [ $1 == untar ] ; then
   for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
   do
       mkdir -p $datadir/$i/raw
       mkdir -p $datadir/$i/batch_logs
       sbatch --output=$datadir/$i/batch_logs/untar.txt $srcdir/untar.scr $datadir/*$i*.tar.gz $datadir/$i
   done
fi

if [ $1 == call ] ; then
    ##for i in antibody mock infected ;
    for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	mkdir -p $datadir/$i/called
	mkdir -p $datadir/$i/call_done

	numdirs=`find $datadir/$i/raw/* -maxdepth 0 -type d | wc -l `
	dummy=1
	maxdir=`expr $numdirs - $dummy`
	echo $maxdir
	sbatch --array=0-$maxdir --output=$datadir/$i/batch_logs/call.txt --job-name=call_$i $srcdir/call_rna.scr $datadir/$i
    done
fi


if [ $1 == fastqs ] ; then
    for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	mkdir -p $datadir/$i/fqs
	cat $datadir/$i/called/*/workspace/pass/*fastq > $datadir/$i/fqs/$i.fq
    done
fi


if [ $1 == align ] ; then
    ml samtools
    for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	mkdir -p $datadir/$i/align

	##sindbis alignment
	minimap2 -a -k14 -uf -t 36 $datadir/refs/sindbis_jane.fasta $datadir/$i/fqs/$i.fq | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/$i/align/$i.sorted.bam -T $datadir/$i/align/reads.tmp -
	samtools index $datadir/$i/align/$i.sorted.bam
	samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.sorted.bam | \
	    samtools sort -@ 36 -o $datadir/$i/align/$i.primary.sorted.bam
	
	##rat spliced alignment
	##minimap2 -a -x splice -uf -k14 -t 36 $datadir/refs/rattus_norvegicus.fa $datadir/$i/fqs/$i.fq  | samtools view -b | samtools sort -o $datadir/$i/align/$i.rat.splicealn.sorted.bam -T $datadir/$i/align/reads.tmp -
	##samtools index $datadir/$i/align/$i.rat.splicealn.sorted.bam
	##samtools view -b -F 0x100 $datadir/$i/align/$i.rat.splicealn.sorted.bam | samtools sort -o $datadir/$i/align/$i.rat.splicealn.primary.sorted.bam
    done
fi


if [ $1 == transfer ] ; then
    for i in Antibody_3dpi Antibody_2dpi Sindbis_2dpi Sindbis_3dpi ;
    do
	scp $datadir/$i/align/$i.primary.sorted.bam smaug:~/Dropbox/Timplab_Data/sindbis/$i/align/$i.primary.sorted.bam
	scp $datadir/$i/align/$i.sorted.bam smaug:~/Dropbox/Timplab_Data/sindbis/$i/align/$i.sorted.bam
    done
fi

	
