#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/mdr
ref=$datadir/dataclean/ref/GRCh38.p14.fa.gz
prefix=200708_mdr_stool16native

if [ $1 == dl_human ] ; then
    mkdir -p $datadir/dataclean
    
    wget -O $ref https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
fi


fq=$datadir/fastqs/$prefix.fq.gz

if [ $1 == alignhuman ] ; then
    minimap2 -t 36 -ax map-ont $ref $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $datadir/dataclean/${prefix}_human.sorted.bam
    samtools index $datadir/dataclean/${prefix}_human.sorted.bam
fi

if [ $1 == unmapped ] ; then
    samtools view -f 4 $datadir/dataclean/${prefix}_human.sorted.bam \
	| awk '{print $1}' \
	      > $datadir/dataclean/nonhuman_readnames.txt
fi

if [ $1 == unmappedfq ] ; then
    seqtk subseq $fq $datadir/dataclean/nonhuman_readnames.txt \
	  > $datadir/dataclean/${prefix}_nohuman.fq
    gzip $datadir/dataclean/${prefix}_nohuman.fq
fi
	
if [ $1 == unmappedf5 ] ; then
    mkdir -p $datadir/dataclean/${prefix}_nohuman

    fast5_subset \
	-i ~/data/mdr/mdr/raw/$prefix/$prefix/20200708_1708_MN19810_FAN26903_ad849ab4/fast5 \
	-s $datadir/dataclean/${prefix}_nohuman \
	-l $datadir/dataclean/nonhuman_readnames.txt \
	-t 36 \
	-r 
fi

hic=181127_hiC_stool
if [ $1 == alignhic ] ; then
    mkdir -p $datadir/dataclean/${hic}_nohuman
    
    ##prep genome
    bwa index $ref

    for i in phase shotgun ;
    do
	bwa mem \
	    -t 36 \
	    $ref \
	    $datadir/illumina/raw/${hic}_${i}_1.fastq.gz \
	    $datadir/illumina/raw/${hic}_${i}_2.fastq.gz |\
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/dataclean/${hic}_${i}.sorted.bam
	samtools index $datadir/dataclean/${hic}_${i}.sorted.bam
    done
fi

    
if [ $1 == unmappedhic ] ; then
    for i in phase shotgun ;
    do
	samtools view -f 4 $datadir/dataclean/${hic}_${i}.sorted.bam \
	    | awk '{print $1}' \
	    | sort \
	    | uniq \
		  > $datadir/dataclean/nonhuman_readnames_$i.txt
    done
fi


if [ $1 == unmappedfqhic ] ; then
    for i in phase shotgun ;
    do
	{
	awk '{print $1"/1"}' $datadir/dataclean/nonhuman_readnames_${i}.txt \
	    > $datadir/dataclean/nonhuman_readnames_${i}_1.txt
	fq1=$datadir/illumina/raw/${hic}_${i}_1.fastq.gz
	seqtk subseq $fq1 $datadir/dataclean/nonhuman_readnames_${i}_1.txt \
	      > $datadir/dataclean/${hic}_nohuman/${hic}_${i}_1_nohuman.fq
	gzip $datadir/dataclean/${hic}_nohuman/${hic}_${i}_1_nohuman.fq

	awk '{print $1"/2"}' $datadir/dataclean/nonhuman_readnames_${i}.txt \
	    > $datadir/dataclean/nonhuman_readnames_${i}_2.txt
	fq2=$datadir/illumina/raw/${hic}_${i}_2.fastq.gz
	seqtk subseq $fq2 $datadir/dataclean/nonhuman_readnames_${i}_2.txt \
	      > $datadir/dataclean/${hic}_nohuman/${hic}_${i}_2_nohuman.fq
	gzip $datadir/dataclean/${hic}_nohuman/${hic}_${i}_2_nohuman.fq 
	} & 
    done
fi
