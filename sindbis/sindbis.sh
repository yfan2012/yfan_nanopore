#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/sindbis
srcdir=~/Code/utils/marcc

if [ $1 == untar ] ; then
   for i in antibody mock infected ;
   do
       mkdir -p $datadir/$i/raw
       mkdir -p $datadir/$i/batch_logs
       sbatch --output=$datadir/$i/batch_logs/untar.txt $srcdir/untar.scr $datadir/$i/*.tar.gz $datadir/$i
   done
fi

if [ $1 == call ] ; then
    for i in antibody mock infected ;
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

if [ $1 == unzip ] ; then
    ##turns out norah already called the runs on gridion
    for i in antibody mock infected ;
    do
	gunzip $datadir/$i/gridfqs/*fq.gz &
    done
fi

if [ $1 == centrifuge ] ; then
    for i in antibody mock infected ;
    do
	dbdir=~/scratch/centrifuge_db
	mkdir -p $datadir/$i/classification

	~/scratch/centrifuge/centrifuge -p 36 -x $dbdir/abv -U $datadir/$i/gridfqs/*.fq -S $datadir/$i/classification/${i}_classification.txt --report-file $datadir/$i/classification/${i}_report.tsv

    done
fi

if [ $1 == kreport ] ; then
    for i in antibody mock infected ;
    do
	dbdir=~/scratch/centrifuge_db
	~/scratch/centrifuge/centrifuge-kreport -x $dbdir/abv $datadir/$i/classification/${i}_classification.txt > $datadir/$i/classification/${i}_kreport.txt
    done
fi

if [ $1 == align_gridfqs ] ; then
    ml samtools
    for i in antibody mock infected ;
    do
	mkdir -p $datadir/$i/align
	minimap2 -a -x map-ont -t 36 $datadir/sindbis_jane.fasta $datadir/$i/gridfqs/*.fq  | samtools view -b | samtools sort -o $datadir/$i/align/$i.gridfq.sorted.bam -T $datadir/$i/align/reads.tmp -
	samtools index $datadir/$i/align/$i.gridfq.sorted.bam
    done
fi


if [ $1 == alignfqs ] ; then
    ml samtools
    ##gatherfq
    for i in antibody mock infected ;
    do
	mkdir -p $datadir/$i/fqs
	cat $datadir/$i/called/*/workspace/pass/*fastq > $datadir/$i/fqs/$i.fq
	minimap2 -a -x map-ont -t 36 $datadir/sindbis_jane.fasta $datadir/$i/fqs/$i.fq  | samtools view -b | samtools sort -o $datadir/$i/align/$i.sorted.bam -T $datadir/$i/align/reads.tmp -
	samtools index $datadir/$i/align/$i.sorted.bam
    done
fi


if [ $1 == cpaligntolocal ] ; then
    ##copy alignments back to local
    for i in antibody mock infected ;
    do
	scp -r $datadir/$i/align smaug:/dilithium/Data/Nanopore/sindbis/$i/
    done
fi


if [ $1 == alignrat ] ; then
    ml samtools
    for i in antibody mock infected ;
    do
	minimap2 -a -x map-ont -t 36 $datadir/rattus_norvegicus.fa $datadir/$i/fqs/$i.fq  | samtools view -b | samtools sort -o $datadir/$i/align/$i.rat.sorted.bam -T $datadir/$i/align/reads.tmp -
	samtools index $datadir/$i/align/$i.rat.sorted.bam
    done
fi

if [ $1 == splicealignrat ] ; then
    ml samtools
    for i in antibody mock infected ;
    do
	minimap2 -a -x splice -uf -k14 -t 36 $datadir/rattus_norvegicus.fa $datadir/$i/fqs/$i.fq  | samtools view -b | samtools sort -o $datadir/$i/align/$i.rat.splicealn.sorted.bam -T $datadir/$i/align/reads.tmp -
	samtools index $datadir/$i/align/$i.rat.sorted.bam
    done
fi
