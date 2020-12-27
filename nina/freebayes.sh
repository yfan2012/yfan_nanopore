#!/bin/bash

##This time, $1 is the results outdir. $2 is the assembly full path. $3 is the prefix. 

mkdir -p $1
mkdir -p $1/index
mkdir -p $1/bam


##copy the raw assembly into the index dir
prefix=$3
cp $2 $1/index/$prefix.fasta

cd ~/software/freebayes/scripts

for i in {1..5} ;
do
    ##build the index and align
    echo building index and aligning for round $i =================================================================
    samtools faidx $1/index/$prefix.fasta
    bowtie2-build $1/index/$prefix.fasta $1/index/$prefix
    bowtie2 -p 54 \
	    -x $1/index/$prefix \
	    -1 $1/*fwd_paired.fastq.gz \
	    -2 $1/*rev_paired.fastq.gz \
	    -r $1/*fwd_unpaired.fastq.gz $1/*rev_unpaired.fastq.gz |\
	samtools view -@ 54 -bS - |\
	samtools sort -@ 54 -o $1/bam/$prefix.sorted.bam
    samtools index $1/bam/$prefix.sorted.bam

    
    ##do the correction
    echo correcting for round $i ================================================================================
    ./freebayes-parallel \
	<(./fasta_generate_regions.py $1/index/$prefix.fasta.fai 100000) 36\
	-f $1/index/$prefix.fasta \
	$1/bam/$prefix.sorted.bam > $1/${prefix}_fb${i}_raw.vcf


    vcffilter -f "AO > RO & AO > 5 & AF > .5" $1/${prefix}_fb${i}_raw.vcf > $1/${prefix}_fb${i}_bwa.vcf
    bgzip -c $1/${prefix}_fb${i}_bwa.vcf > $1/${prefix}_fb${i}_bwa.vcf.gz
    tabix -p vcf $1/${prefix}_fb${i}_bwa.vcf.gz
    bcftools consensus $1/${prefix}_fb${i}_bwa.vcf.gz < $1/index/$prefix.fasta > $1/${prefix}_fb${i}_bwa.fasta

    ##newly corrected genome replaces the old genome in the index dir
    echo moving old
    mv $1/bam/$prefix.sorted.bam $1/bam/$prefix.$i.sorted.bam
    mv $1/bam/$prefix.sorted.bam.bai $1/bam/$prefix.$i.sorted.bam.bai
    rm $1/index/*
    echo copying $i to empty index folder
    cp $1/${prefix}_fb${i}_bwa.fasta $1/index/$prefix.fasta

done


