#!/bin/bash

##This time, $1 is the results outdir. $2 is the assembly full path. $3 is the prefix. 

mkdir -p $1
mkdir -p $1/index
mkdir -p $1/bam


##copy the raw assembly into the index dir
prefix=$3
cp $2 $1/index/$prefix.fasta


for i in {1..3} ;
do
    ##build the index and align
    echo building index and aligning for round $i =================================================================
    samtools faidx $1/index/$prefix.fasta
    bwa index $1/index/$prefix.fasta
    bwa mem -t 36 $1/index/$prefix.fasta $1/*fwd_paired.fq.gz $1/*rev_paired.fq.gz |\
	samtools view -@ 36 -bS - |\
	samtools sort -@ 36 -o $1/bam/$prefix.sorted.bam
    samtools index $1/bam/$prefix.sorted.bam

     
    ##do the correction
    echo correcting for round $i ================================================================================
    ./freebayes-parallel \
	<(./fasta_generate_regions.py $1/index/$prefix.fasta.fai 100000) 36\
	-f $1/index/$prefix.fasta \
	$1/bam/$prefix.sorted.bam > $1/nivar_fb${i}_bwa_raw.vcf


    vcffilter -f "AO > 3"  $1/nivar_fb${i}_bwa_raw.vcf > $1/nivar_fb${i}_bwa.vcf
    bgzip -c $1/nivar_fb${i}_bwa.vcf > $1/nivar_fb${i}_bwa.vcf.gz
    tabix -p vcf $1/nivar_fb${i}_bwa.vcf.gz
    bcftools consensus $1/nivar_fb${i}_bwa.vcf.gz < $1/index/$prefix.fasta > $1/nivar_fb${i}_bwa.fasta

    ##see how many changes
    dnadiff $1/index/$prefix.fasta $1/nivar_fb${i}_bwa.fasta -d $i/nivar_fb${i}.delta
    
    ##newly corrected genome replaces the old genome in the index dir
    echo clearing old
    rm $1/index/*
    rm $1/bam/*
    echo copying $i to empty index folder
    cp $1/nivar_fb${i}_bwa.fasta $1/index/$prefix.fasta

    
done


