#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans
dbxdir=~/Dropbox/timplab_data/prolificans

if [ $1 == mummer ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/mummer

	nucmer \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/genomes/$i.flye.fasta \
	    $datadir/$i/genomes/$i.canu.fasta
	
	mummerplot \
	    --filter --fat --postscript \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/mummer/$i.delta
	
	mummerplot \
	    --filter --fat --png \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/mummer/$i.delta

	dnadiff \
	    -p $datadir/$i/mummer/$i \
	    $datadir/$i/genomes/$i.flye.fasta \
	    $datadir/$i/genomes/$i.canu.fasta
    done
fi

if [ $1 == asmstats ] ; then

    mkdir -p $dbxdir/qc
    touch $dbxdir/qc/asmstats.csv
    
    for i in st31 st90853 st5317 ;
    do
	python ~/Code/utils/qc/asm_assess.py \
	       -i $datadir/$i/genomes/$i.canu.fasta \
	       -p $i.canu >> $dbxdir/qc/asmstats.csv
	python ~/Code/utils/qc/asm_assess.py \
	       -i $datadir/$i/genomes/$i.flye.fasta \
	       -p $i.flye >> $dbxdir/qc/asmstats.csv
	python ~/Code/utils/qc/asm_assess.py \
	       -i $datadir/$i/ragtag/cf/$i.ragtag_cf.scaffolds.fasta \
	       -p $i.ragtag_cf >> $dbxdir/qc/asmstats.csv
	python ~/Code/utils/qc/asm_assess.py \
	       -i $datadir/$i/ragtag/fc/$i.ragtag_fc.scaffolds.fasta \
	       -p $i.ragtag_fc >> $dbxdir/qc/asmstats.csv
    done
fi


if [ $1 == runstats ] ; then

    touch $dbxdir/qc/runstats.csv
    touch $dbxdir/qc/runstats_long.csv
    
    for i in st31 st90853 st5317 ;
    do
	bash ~/Code/utils/qc/basic_run_assess.sh \
	     $datadir/$i/reads/ont/${i}.fastq.gz >> $dbxdir/qc/runstats.csv
	bash ~/Code/utils/qc/basic_run_assess.sh \
	     $datadir/$i/reads/ont/${i}_long.fastq.gz >> $dbxdir/qc/runstats_long.csv
    done
fi



if [ $1 == align ] ; then

    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/align
	fq=$datadir/$i/reads/ont/${i}_long.fastq.gz

	for asm in $datadir/$i/genomes/*fasta ;
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
	
	for asm in $datadir/$i/genomes/*fasta ;
	do
	    prefix=`basename $asm .fasta`
	    
	    bowtie2-build -q $asm $datadir/$i/genomes/$prefix
	    bowtie2 -p 54 \
		    -x $datadir/$i/genomes/$prefix \
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
	
	for align in $datadir/$i/align/*.sorted.bam ;
	do
	    prefix=`basename $align .sorted.bam`
	    refprefix=`echo $prefix | cut -d . -f 1,2`
	    echo $refprefix
	    ref=$datadir/$i/genomes/$refprefix.fasta
	    samtools faidx $ref
	    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $ref.fai > $ref.bed
	    bedtools coverage -d \
		     -a $ref.bed \
		     -b $align > $datadir/$i/cov/$prefix.cov
	done
    done
fi

if [ $1 == find_breaks ] ; then
    Rscript ./qc.R
fi

if [ $1 == breaktigs ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/genomes_covfilt
	for genome in $datadir/$i/genomes/*fasta ;
	do
	    prefix=`basename $genome .fasta`
	    grep $prefix $dbxdir/zero_cov_telofilt.csv > $dbxdir/zero_cov_telofilt_tmp.csv

	    python ~/Code/yfan_nanopore/nina/breaks.py \
		   -a $genome \
		   -n $datadir/$i/align/$prefix.sorted.bam \
		   -i $datadir/$i/align/$prefix.illumina.sorted.bam \
		   -r $dbxdir/zero_cov_telofilt_tmp.csv \
		   -o $datadir/$i/genomes_covfilt/$prefix.covfilt.fasta
	done
    done
fi
	    

##giving up on denovo mito discovery: just grab from laurent's circ asm
ref=$datadir/ref/LProlificans_v1.0.fa
if [ $1 == grabmito ] ; then
    samtools faidx $ref
    bytestart=`grep mitoscaff1 $ref.fai | awk '{print $3}'`
    length=`grep mitoscaff1 $ref.fai | awk '{print $2}'`
    linelen=`grep mitoscaff1 $ref.fai | awk '{print $4}'`
    numlines=$(( length/linelen ))
    finlines=$(( $numlines + 2 ))
    grep mitoscaff1 $ref | tr '\n' ' ' >> $datadir/ref/LProlificans_mito_v1.0.fa
    tail -c +$bytestart $ref | head -n $finlines >> $datadir/ref/LProlificans_mito_v1.0.fa
fi
  

if [ $1 == mummer_mito ] ; then
    mito=$datadir/ref/LProlificans_mito_v1.0.fa
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/mummer_mito

	for genome in $datadir/$i/genomes_covfilt/*fasta ;
	do
	    prefix=`basename $genome .covfilt.fasta`

	    nucmer \
		-p $datadir/$i/mummer_mito/$prefix \
		$genome \
		$mito
	    
	    mummerplot \
		--filter --fat --postscript \
		-p $datadir/$i/mummer_mito/$prefix \
		$datadir/$i/mummer_mito/$prefix.delta
	    
	    mummerplot \
		--filter --fat --png \
		-p $datadir/$i/mummer_mito/$prefix \
		$datadir/$i/mummer_mito/$prefix.delta
	    
	    dnadiff \
		-p $datadir/$i/mummer_mito/$prefix \
		$genome \
		$mito
	done
    done
fi


if [ $1 == mito_trim ] ; then
    Rscript find_mito.R
fi



if [ $1 == final ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/final
    done

    cp $datadir/st90853/genomes_final/st90853.ragtag_fc.final.fasta $datadir/st90853/final/st90853.final.fasta
    cp $datadir/st31/genomes_final/st31.ragtag_fc.final.fasta $datadir/st31/final/st31.final.fasta
    cp $datadir/st5317/genomes_final/st5317.ragtag_cf.final.fasta $datadir/st5317/final/st5317.final.fasta

    for i in st31 st90853 st5317 ;
    do
	asm=$datadir/$i/final/$i.final.fasta
	sed -i -e 's/000000//g' $asm
	sed -i -e "s/st31.canu/$i/g" $asm
    done
fi


if [ $1 == final_busco ] ; then
    mkdir -p $datadir/st90853/final/busco
    cp -r $datadir/st90853/busco/st90853.ragtag_fc/st90853.ragtag_fc/* $datadir/st90853/final/busco/

    mkdir -p $datadir/st31/final/busco
    cp -r $datadir/st31/busco/st31.ragtag_fc/st31.ragtag_fc/* $datadir/st31/final/busco/

    mkdir -p $datadir/st5317/final/busco
    cp -r $datadir/st5317/busco/st5317.ragtag_cf/st5317.ragtag_cf/* $datadir/st5317/final/busco/
fi


    
