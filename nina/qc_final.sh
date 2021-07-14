#!/bin/bash

##cov based breaking, mito trimming, etc is done

datadir=/pym/Data/Nanopore/projects/prolificans
dbxdir=~/Dropbox/timplab_data/prolificans

if [ $1 == trimmed_asmstats ] ; then
    mkdir -p $dbxdir/qc
    touch $dbxdir/qc/asmstats_trimmed.csv
    
    for i in st31 st90853 st5317 ;
    do
	for asm in $datadir/$i/genomes_mitotrim/*fasta ;
	do
	    prefix=`basename $asm .fasta`
	    python ~/Code/utils/qc/asm_assess.py \
		   -i $asm \
		   -p $prefix >> $dbxdir/qc/asmstats_trimmed.csv
	done
    done
fi

if [ $1 == align ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/align
	fq=$datadir/$i/reads/ont/${i}_long.fastq.gz

	for asm in $datadir/$i/genomes_mitotrim/*fasta ;
	do
	    prefix=`basename $asm .fasta`
	    minimap2 -t 54 -ax map-ont $asm $fq |\
		samtools view -@ 54 -b |\
		samtools sort -@ 54 -o $datadir/$i/align/$prefix.sorted.bam
	    samtools index $datadir/$i/align/$prefix.sorted.bam
	done
    done
fi

if [ $1 == align_illumina ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/align
	r1=$datadir/$i/reads/illumina/${i}_fwd_paired.fq.gz
	r2=$datadir/$i/reads/illumina/${i}_rev_paired.fq.gz
	
	for asm in $datadir/$i/genomes_mitotrim/*fasta ;
	do
	    prefix=`basename $asm .fasta`
	    
	    bowtie2-build -q $asm $datadir/$i/genomes_mitotrim/$prefix
	    bowtie2 -p 54 \
		    -x $datadir/$i/genomes_mitotrim/$prefix \
		    -1 $r1 \
		    -2 $r2 |\
		samtools view -@ 54 -b |\
		samtools sort -@ 54 -o $datadir/$i/align/$prefix.illumina.sorted.bam
	    samtools index $datadir/$i/align/$prefix.illumina.sorted.bam
	done
    done
fi

if [ $1 == coverage ] ;then
    for i in st31 st90853 st5317 ;
    do
        mkdir -p $datadir/$i/cov

        for align in $datadir/$i/align/*mitotrim*.sorted.bam ;
        do
            prefix=`basename $align .sorted.bam`
            refprefix=`echo $prefix | cut -d . -f 1,2,3`
            echo $refprefix
            ref=$datadir/$i/genomes_mitotrim/$refprefix.fasta
            samtools faidx $ref
            awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $ref.fai > $ref.bed
            bedtools coverage -d \
                     -a $ref.bed \
                     -b $align > $datadir/$i/cov/$prefix.cov
        done
    done
fi

if [ $1 == make_final ] ; then
    Rscript qc_final.R
fi


if [ $1 == final_stats ] ; then
    mkdir -p $dbxdir/qc
    touch $dbxdir/qc/asmstats_final.csv
    
    for i in st31 st90853 st5317 ;
    do
	for asm in $datadir/$i/genomes_final/*fasta ;
	do
	    prefix=`basename $asm .fasta`
	    python ~/Code/utils/qc/asm_assess.py \
		   -i $asm \
		   -p $prefix >> $dbxdir/qc/asmstats_final.csv
	done
    done
fi

if [ $1 == busco ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/busco

	for asm in $datadir/$i/genomes_final/*fasta ;
	do
	    prefix=`basename $asm .final.fasta`
	    mkdir -p $datadir/$i/busco/$prefix
	    busco \
		-m genome \
		-l sordariomycetes_odb10 \
		-i $asm \
		-o $prefix \
		--out_path $datadir/$i/busco/$prefix \
		-c 36 \
		-f
	done
	
    done
fi

if [ $1 == ref_busco ] ; then
    busco \
	-m genome \
	-l sordariomycetes_odb10 \
	-i $datadir/ref/LProlificans_v1.0.fa \
	-o LProlificans_v1.0 \
	--out_path $datadir/ref/busco_LProlificans_v1.0 \
	-c 36 \
	-f
fi

    
		     

    
