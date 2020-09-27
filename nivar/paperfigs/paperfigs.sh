#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=/uru/Data/Nanopore/projects/nivar/paperfigs
fq=/uru/Data/Nanopore/projects/nivar/r9/r9_3kb.fq
dbxdir=~/Dropbox/yfan/nivar/paperfigs

if [ $1 == runqc ] ; then
    Rscript ~/Code/utils/qc/run_summary.R \
	    -i $rawdir/r9/run1_called/sequencing_summary.txt \
	    -o $dbxdir/run1_run_summary.pdf \
	    -p Run1
    Rscript ~/Code/utils/qc/run_summary.R \
	    -i $rawdir/r9/run2_called/sequencing_summary.txt \
	    -o $dbxdir/run2_run_summary.pdf \
	    -p Run2
    Rscript ~/Code/utils/qc/run_summary.R \
	    -i $rawdir/r10/called/sequencing_summary.txt \
	    -o $dbxdir/r10_summary.pdf \
	    -p r10
fi

if [ $1 == assemble ] ; then
    mkdir $datadir/assembly
    
    canu \
	-p nivar \
	-d $datadir/assembly \
	genomeSize=11.6m \
	-nanopore-raw $fq
fi

if [ $1 == medaka ] ; then
    mkdir -p $datadir/medaka
    
    medaka_consensus \
	-i $fq \
	-d $datadir/assembly/nivar.contigs.fasta \
	-o $datadir/medaka \
	-t 36 \
	-m r941_min_high_g303
fi


if [ $1 == polish ] ; then
    mkdir -p $datadir/freebayes
    cp $rawdir/illumina/gDNA_trimmed/*_paired.fq.gz $datadir/freebayes
    
    bash ./freebayes_bwa.sh $datadir/freebayes $datadir/assembly/nivar.contigs.fasta nivar
fi

if [ $1 == medusa ] ; then
    mkdir -p $datadir/medusa
    cp -r ~/software/medusa/medusa_scripts ./

    java -jar ~/software/medusa/medusa.jar \
	 -f $rawdir/reference/medusa_fungi \
	 -i $datadir/freebayes/nivar_fb3_bwa.fasta \
	 -v \
	 -o $datadir/medusa/nivar.scaffold.fasta
fi

ref=$rawdir/reference/candida_nivariensis.fa

if [ $1 == busco ] ; then
    python ~/software/busco/scripts/run_BUSCO.py -f \
	   -i $datadir/medusa/nivar.scaffold.fasta \
	   -l ~/software/busco/lineages/fungi_odb9 \
	   -sp candida_albicans \
	   -o nivar \
	   -m genome
    python ~/software/busco/scripts/run_BUSCO.py -f \
	   -i $ref \
	   -l ~/software/busco/lineages/fungi_odb9 \
	   -sp candida_albicans \
	   -o ref \
	   -m genome

    mv ./run_* $datadir/busco/
fi
