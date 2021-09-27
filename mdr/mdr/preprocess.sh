#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/mdr
prefix=200708_mdr_stool16native

ref=$datadir/ref/mdr_refs.fa
asm=$datadir/flye/$prefix/$prefix.assembly.fasta
asmpolished=$datadir/medaka/consensus.fasta

if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs

    cat $datadir/called/pass/*fastq.gz > $datadir/fastqs/$prefix.fq.gz
fi

if [ $1 == flye ] ; then
    mkdir -p $datadir/flye
    mkdir -p $datadir/flye/$prefix

    flye \
	--nano-raw $datadir/fastqs/$prefix.fq.gz \
	-o $datadir/flye/$prefix \
	-t 36 \
	-g 100m \
	--meta
    
    mv $datadir/flye/$prefix/assembly.fasta $datadir/flye/$prefix/$prefix.assembly.fasta
fi

fq=$datadir/fastqs/200708_mdr_stool16native.fq.gz
if [ $1 == racon_align ] ; then
    mkdir -p $datadir/racon

    minimap2 -t 36 -ax map-ont $datadir/flye/$prefix/$prefix.assembly.fasta $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $datadir/racon/${prefix}_asm.sorted.bam
    samtools index $datadir/racon/${prefix}_asm.sorted.bam
fi

if [ $1 == racon ] ; then
    mkdir -p $datadir/racon

    samtools view -@ 36 $datadir/racon/${prefix}_asm.sorted.bam > $datadir/racon/${prefix}_asm.sorted.sam
    racon -m 8 -x -6 -g -8 -w 500 -t 54\
	  $fq \
	  $datadir/racon/${prefix}_asm.sorted.sam \
	  $datadir/flye/$prefix/$prefix.assembly.fasta > $datadir/racon/$prefix.racon.contigs.fasta
fi

if [ $1 == medaka ] ; then
    mkdir -p $datadir/medaka

    medaka_consensus \
	-i $fq \
	-d $datadir/racon/$prefix.racon.contigs.fasta \
	-o $datadir/medaka \
	-t 54 \
	-m r941_min_high_g360
fi

if [ $1 == amr ] ; then
    mkdir -p $datadir/amr

    for i in ecoh card ncbi resfinder plasmidfinder vfdb ecoli_vf megares argannot ;
    do
	abricate \
	    --threads 36 \
	    --db $i \
	    $asmpolished > $datadir/amr/$prefix.$i.tsv
    done
fi


if [ $1 == blastreads ] ; then
    mkdir -p $datadir/blast_reads

    seqtk seq -a $datadir/fastqs/$prefix.fq.gz > $datadir/fastqs/$prefix.fasta
    
    blastn \
	-num_threads 36 \
	-query $datadir/fastqs/$prefix.fasta \
	-db /atium/Data/ref/ncbi/nt \
	-outfmt 7 \
	-out $datadir/blast_reads/$prefix.reads.tsv
fi

if [ $1 == blastflye ] ; then
    mkdir -p $datadir/blast_contigs

    blastn \
	-num_threads 36 \
	-query $datadir/flye/$prefix/$prefix.assembly.fasta \
	-db /atium/Data/ref/ncbi/nt \
	-outfmt 7 \
	-out $datadir/blast_contigs/$prefix.assembly.tsv
fi

if [ $1 == blastpolished ] ; then
    mkdir -p $datadir/blast_contigs
    blastn \
	-num_threads 36 \
	-query $asmpolished \
	-db /atium/Data/ref/bacteria/blastdb/all_bacteria_refs \
	-outfmt 7 \
	-out $datadir/blast_contigs/${prefix}_polished.assembly.tsv
fi

if [ $1 == blastpolished_ntdb ] ; then
    mkdir -p $datadir/blast_contigs
    blastn \
	-num_threads 36 \
	-query $asmpolished \
	-db /atium/Data/ref/ncbi/nt \
	-outfmt 7 \
	-out $datadir/blast_contigs/${prefix}_ntdb_polished.assembly.tsv
fi    

dbdir=/mithril/Data/Nanopore/ref/kraken2

if [ $1 == kraken ] ; then
    mkdir -p $datadir/kraken

    fq=$datadir/fastqs/$prefix.fq.gz

    kraken2 \
	--db $dbdir/standard_ont \
	--threads 36 \
	--classified-out $datadir/kraken/$prefix.class.txt \
	--unclassified-out $datadir/kraken/$prefix.unclass.txt \
	--output $datadir/kraken/$prefix.out.txt \
	--report $datadir/kraken/$prefix.report.txt \
	--use-names \
	$fq
fi

if [ $1 == top40 ] ; then
    awk '$4 == "S" {print $0}' $datadir/kraken/$prefix.report.txt | \
	sort -r -k2 -n | \
	head -n 40 > $datadir/kraken/$prefix.report.top40.txt
    awk '$4 == "P" {print $0}' $datadir/kraken/$prefix.report.txt | \
	sort -r -k2 -n | \
	head -n 40 > $datadir/kraken/$prefix.report.top40_phylum.txt
fi


if [ $1 == cattigs ] ; then
    catdir=/atium/Data/ref/CATdb
    mkdir -p $datadir/contig_id/CAT

    CAT contigs \
	-c $datadir/medaka/consensus.fasta \
	-d $catdir/db \
	-t $catdir/tax \
	-o $datadir/contig_id/CAT/$prefix.CAT \
	--force \
	--sensitive
    
fi
