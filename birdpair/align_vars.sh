#!/bin/bash

##annotate variants based on alignment

datadir=/kyber/Data/Nanopore/projects/birdpair

if [ $1 == align ] ; then
    mkdir -p $datadir/align

    for i in $datadir/ref/*fa ;
    do
	refname=`basename $i .fa`
	bwa index $i

	##don't include unmapped, non-primary, failqual, pcr dup, supplementary
	bwa mem -t 36 $i $datadir/trimmed/patient_fwd_paired.fq.gz $datadir/trimmed/patient_rev_paired.fq.gz |
	    samtools view -@ 36 -b -q 30 -F 3844 |
	    samtools sort -@ 36 -o $datadir/align/patient_$refname.sorted.bam
	samtools index $datadir/align/patient_$refname.sorted.bam

	bwa mem -t 36 $i $datadir/trimmed/bird_fwd_paired.fq.gz $datadir/trimmed/bird_rev_paired.fq.gz |
	    samtools view -@ 36 -b -q 30 -F 3844 |
	    samtools sort -@ 36 -o $datadir/align/bird_$refname.sorted.bam
	samtools index $datadir/align/bird_$refname.sorted.bam

    done
fi

    
if [ $1 == freebayes ] ; then
    mkdir -p $datadir/vars
    
    
    cd ~/software/freebayes/scripts
    for i in $datadir/align/*sorted.bam ;
    do
	ref=`basename $i .sorted.bam | cut -d _ -f 2,3`
	samtools faidx $datadir/ref/$ref.fa

	prefix=`basename $i .sorted.bam`
	./freebayes-parallel \
	    <(./fasta_generate_regions.py $datadir/ref/$ref.fa.fai 100000) 36\
	    -f $datadir/ref/$ref.fa \
	    $datadir/align/$prefix.sorted.bam > $datadir/vars/${prefix}_raw.vcf

	vcffilter -f "AO > 3 & RO < 1"  $datadir/vars/${prefix}_raw.vcf > $datadir/vars/${prefix}.vcf
    done
fi
	
if [ $1 == annot ] ; then


    for i in $datadir/vars/*vcf ;
    do
	prefix=`basename $i .vcf`
	ref=`echo $prefix | cut -d _ -f 2,3`
	python ~/Code/utils/find_snps.py -v $i -g $datadir/ref/$ref.gff -o $datadir/vars/$prefix.csv
    done
fi
	     
    
