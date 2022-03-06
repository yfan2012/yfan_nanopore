#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/barnyard/strains/bisulfite
refdir=/mithril/Data/Nanopore/projects/methbin/barnyard/ref/bisulf
ref=$refdir/allrefs.fa

prefix=220131_mdr_barnyard_
samps='st3294 st3689'

if [ $1 == prep_genome ] ; then
    bismark_genome_preparation \
	--parallel 4 \
	$refdir
fi

if [ $1 == bismark ] ; then
    mkdir -p $datadir/bismark
    for i in $samps ;
    do
	bismark \
	    --genome $refdir \
	    --parallel 12 \
	    -1 $datadir/raw/*${i}*R1*.fastq.gz \
	    -2 $datadir/raw/*${i}*R2*.fastq.gz \
	    -o $datadir/bismark/$i
    done
fi


if [ $1 == extractor ] ; then
    for i in $samps ;
    do
	alignname=`basename $datadir/bismark/$i/*bam .bam`
	samtools view -@ 36 $datadir/bismark/$i/$alignname.bam > $datadir/bismark/$i/$alignname.sam
	bismark_methylation_extractor \
	    --comprehensive \
	    --merge_non_CpG \
	    --parallel 12 \
	    --cytosine_report \
	    --CX \
	    --genome_folder $refdir \
	    -o $datadir/bismark/$i \
	    $datadir/bismark/$i/$alignname.bam
    done
fi


barcodes=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes2.txt
if [ $1 == count ] ; then
    mkdir -p $datadir/motifcalls
    for i in $samps ;
    do
	for seqname in PRW62 NC_007795.1 CP017100.1 ;
	do
	    python3 ~/Code/yfan_meth/utils/bismark_motif_finder.py \
		    -c $datadir/bismark/$i/*$i*.CX_report.txt \
		    -r $ref \
		    -b $barcodes \
		    -p .7 \
		    -l 8 \
		    -m 3 \
		    -o $datadir/motifcalls/$i.$seqname.seqs.fasta \
		    -s $seqname \
		    >> $datadir/motifcalls/$i.clinbc2.csv
	done
    done
fi

	       
