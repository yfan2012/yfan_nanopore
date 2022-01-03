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

    
if [ $1 == getchrlist ] ; then
    ##get list of bacterial contigs from merged reference
    ##just going to manually add the sample labels

    grep '>' $ref \
	| awk '{print substr($1,2)}' \
	| awk ' { if (substr($1,1,1)!="t" ) print $0}' \
	      > ./chrlist.txt
fi
    

if [ $1 == test_bisulf ] ; then

    while read p; do
	label=`echo $p | cut -d ' ' -f 2`
	chrom=`echo $p | cut -d ' ' -f 1`
	python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
	       -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
	       -r $ref \
	       -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
	       -p .5 \
	       -l 8 \
	       -m 10 \
	       -s $chrom
    done <chrlist_test.txt
fi

if [ $1 == bisulf ] ; then
    mkdir -p $datadir/motifcalls
    touch $datadir/motifcalls/50_barcodes.csv
    while read p; do
        label=`echo $p | cut -d ' ' -f 2`
        chrom=`echo $p | cut -d ' ' -f 1`
        python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
               -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
               -p .5 \
               -l 8 \
               -m 10 \
               -s $chrom >> $datadir/motifcalls/50_barcodes.csv
    done <chrlist_withlabels.txt
fi


if [ $1 == pseudo ] ; then
    label=paeruginosa
    chrom=Pseudomonas_aeruginosa_complete_genome
    python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
           -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
           -p .5 \
           -l 8 \
           -m 10 \
	   -o $datadir/motifcalls/$label.seqs.fasta \
           -s $chrom >> $datadir/motifcalls/50_barcodes.csv
fi

if [ $1 == ecoli_investigation ] ; then
    label=ecoli
    chrom=Escherichia_coli_chromosome
    python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
           -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
           -r $ref \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
           -p .5 \
           -l 8 \
           -m 10 \
	   -o $datadir/motifcalls/$label.seqs.fasta \
           -s $chrom >> $datadir/motifcalls/50_barcodes.csv
fi


if [ $1 == bisulf_pseudo ] ; then
    ##rerun with pseudomonas motif
    mkdir -p $datadir/motifcalls
    touch $datadir/motifcalls/51_barcodes.csv
    while read p; do
        label=`echo $p | cut -d ' ' -f 2`
        chrom=`echo $p | cut -d ' ' -f 1`
        python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
               -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes51.txt \
               -p .5 \
               -l 8 \
               -m 10 \
               -s $chrom >> $datadir/motifcalls/51_barcodes.csv
    done <chrlist_withlabels.txt
fi

if [ $1 == countunk ] ; then
    ##rerun with pseudomonas motif
    mkdir -p $datadir/motifcalls
    touch $datadir/motifcalls/51_barcodes.csv
    while read p; do
        label=`echo $p | cut -d ' ' -f 2`
        chrom=`echo $p | cut -d ' ' -f 1`
        python ~/Code/yfan_meth/utils/bismark_motif_finder.py \
               -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
               -r $ref \
               -b ~/Code/yfan_nanopore/mdr/rebase/barcodes51.txt \
               -p .5 \
               -l 8 \
               -m 10 \
	       -u \
               -s $chrom >> $datadir/motifcalls/51_barcodes_unk.csv
    done <chrlist_withlabels.txt
fi

if [ $1 == account ] ; then 
    while read p; do
	label=`echo $p | awk '{print $1}'`
	motifs=`echo $p | awk '{print $2}'`
	methpos=`echo $p | awk '{print $3}'`
	
	python ~/Code/yfan_meth/utils/bismark_motif_confirm.py \
	       -c $datadir/bismark/$label/${label}_1_bismark_bt2_pe.CX_report.txt \
	       -r $ref \
	       -m $motifs \
	       -p $methpos
    done <bismark_account.txt
fi
