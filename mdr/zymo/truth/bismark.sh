#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite
refdir=/uru/Data/Nanopore/projects/read_class/zymo/ref
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == prep_genome ] ; then
    ##do genome prep step
    ##most settings should be okay with this application
    bismark_genome_preparation \
	--parallel 4 \
	$refdir
fi

if [ $1 == bismark ] ; then
    mkdir -p $datadir/bismark
    for i in bsubtilis ecoli efaecalis lmonocytogenes nasa_Ecoli_K12 paeruginosa saureus senterica ;
    do

	bismark \
	    --genome $refdir \
	    --parallel 8 \
	    -1 $datadir/raw/fastqs/${i}_1.fastq \
	    -2 $datadir/raw/fastqs/${i}_2.fastq \
	    -o $datadir/bismark/$i
    done
fi


if [ $1 == extractor ] ; then
    for i in bsubtilis ecoli efaecalis lmonocytogenes nasa_Ecoli_K12 paeruginosa saureus senterica ;
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

    
     
