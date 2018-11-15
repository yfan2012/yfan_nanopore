#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/sindbis
srcdir=~/Code/utils/marcc

if [ $1 == backup ] ; then
    for i in antibody mock infected;
    do
	aws s3 cp $datadir/$i/*.tar.gz s3://yfan-seqruns/griffin/
    done
fi

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

	~/scratch/centrifuge/centrifuge -p 36 -x $dbdir/abvm -U $datadir/$i/fqs/*.fq -S $datadir/$i/classification/${i}_classification.txt --report-file $datadir/$i/classification/${i}_report.tsv
    done
fi

if [ $1 == kreport ] ; then
    for i in antibody mock infected ;
    do
	dbdir=~/scratch/centrifuge_db
	~/scratch/centrifuge/centrifuge-kreport -x $dbdir/abvm $datadir/$i/classification/${i}_classification.txt > $datadir/$i/classification/${i}_kreport.txt
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

if [ $1 == cpclasstolocal ] ; then
    ##copy alignments back to local
    for i in antibody mock infected ;
    do
	scp -r $datadir/$i/classification smaug:/dilithium/Data/Nanopore/sindbis/$i/
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


if [ $1 == transcript_alignrat ] ; then
    ml samtools
    for i in antibody mock infected ;
    do
	minimap2 -a -x splice -uf -k14 -t 24 $datadir/rattus_norvegicus.rna.fa $datadir/$i/fqs/$i.fq  | samtools view -b | samtools sort -o $datadir/$i/align/$i.rat.transcriptaln.sorted.bam -T $datadir/$i/align/reads.tmp -
	samtools index $datadir/$i/align/$i.rat.transcriptaln.sorted.bam
    done
fi


if [ $1 == get_primary_only ] ; then
    ml samtools
    for i in antibody mock infected ;
	do
	(samtools view -b -F 0x100 $datadir/$i/align/$i.rat.splicealn.sorted.bam | samtools sort -o $datadir/$i/align/$i.rat.splicealn.primary.sorted.bam
	samtools view -b -F 0x100 $datadir/$i/align/$i.rat.transcriptaln.sorted.bam | samtools sort -o $datadir/$i/align/$i.rat.transcriptaln.primary.sorted.bam
	samtools view -b -F 0x100 $datadir/$i/align/$i.sorted.bam | samtools sort -o $datadir/$i/align/$i.primary.sorted.bam) &
    done
fi


if [ $1 == count_rat_transcripts ] ; then
    for i in antibody mock infected ;
    do
	mkdir -p $datadir/$i/rat_transcripts
	python rat_expression.py -i $datadir/$i/align/$i.rat.transcriptaln.primary.sorted.bam -r $datadir/rattus_norvegicus.rna.fa -o $datadir/$i/rat_transcripts/$i.rat_transcript_counts.csv
    done
fi

