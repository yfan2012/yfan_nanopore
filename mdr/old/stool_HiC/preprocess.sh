#!/bin/bash

datadir=/kyber/Data/NGS/projects/181127_hiC_stool

if [ $1 == bmtagger_bitmask ] ; then
   bmtool -d $datadir/bmtagger/GRCH38.fa -o $datadir/bmtagger/GRCH38.bitmask -A 0 -w 18
fi

if [ $1 == bmtagger_srprism ] ; then
    srprism mkindex -i $datadir/bmtagger/GRCH38.fa -o $datadir/bmtagger/GRCH38.srprism -M 7168
fi

if [ $1 == bmtagger_blast ] ; then
    makeblastdb -in $datadir/bmtagger/GRCH38.fa -dbtype nucl
fi

if [ $1 == bmtagger_pefq ] ; then
    for i in $datadir/*R1* ;
    do
	(prefix=` basename $i _R1.fastq`
	bmfilter -b $datadir/bmtagger/GRCH38.bitmask -T $datadir/tmp -q 1 -1 $i -2 $datadir/${prefix}_R2.fastq -o $datadir/$prefix.humanfilt.fastq) &
    done
fi


if [ $1 == bmtagger_script ] ; then
    for i in $datadir/*R1* ;
    do
	prefix=` basename $i _R1.fastq`
	mkdir -p $datadir/tmp
	bmtagger.sh -X -b $datadir/bmtagger/GRCH38.bitmask -x $datadir/bmtagger/GRCH38.srprism -T $datadir/tmp -q 1 -1 $i -2 $datadir/${prefix}_R2.fastq -o $datadir/$prefix.humanfilt.fastq
    done
fi

if [ $1 == align_check ] ; then
    mkdir -p $datadir/align
    ref=$datadir/bmtagger/GRCH38.fa    
    bowtie2-build $ref $datadir/align/GRCH38
    for i in 181127_hiC_stool_shotgun 181127_hiC_stool_phase ;
    do
	bowtie2 --threads 12 -x $datadir/align/GRCH38 -1 $datadir/$i.humanfilt.fastq_1.fastq -2 $datadir/$i.humanfilt.fastq_2.fastq | samtools view -bS | samtools sort -o $datadir/align/$i.humanfilt.sorted.bam
	bowtie2 --threads 12 -x $datadir/align/GRCH38 -1 $datadir/${i}_R1.fastq -2 $datadir/${i}_R2.fastq | samtools view -bS | samtools sort -o $datadir/align/$i.sorted.bam
    done
fi

	

