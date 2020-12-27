#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == grabont ] ; then
    mkdir -p $datadir/st5317/reads/ont

    fastq-dump \
	-O $datadir/st5317/reads/ont \
	--gzip \
	SRR5812844

    mv $datadir/st5317/reads/ont/SRR5812844.fastq.gz $datadir/st5317/reads/ont/st5317.fastq.gz
fi
	
if [ $1 == grabillumina ] ; then
    mkdir -p $datadir/st5317/reads/illumina
    
    fastq-dump \
	-O $datadir/st5317/reads/illumina \
	--gzip \
	--split-files \
	SRR5812813
    
    mv $datadir/st5317/reads/illumina/SRR5812813_1.fastq.gz $datadir/st5317/reads/illumina/st5317_1.fastq.gz
    mv $datadir/st5317/reads/illumina/SRR5812813_2.fastq.gz $datadir/st5317/reads/illumina/st5317_2.fastq.gz
fi

if [ $1 == basecall ] ; then
    ##winry docker
    datadir=/media
    for i in 20180827_2011_nina_fungus2 20181108_0055_181107_nina_v2 ;
    do
	guppy_basecaller \
	    --recursive \
	    --compress_fastq \
	    --kit SQK-LSK109 \
	    --barcode_kits EXP-NBD103 \
	    --trim_barcodes \
	    --flowcell FLO-MIN106 \
	    -i $datadir/$i/fast5 \
	    -s $datadir/$i/called \
	    -x 'cuda:0'
    done
fi

if [ $1 == fqs ] ; then
    mkdir -p $datadir/st90853/reads/ont
    mkdir -p $datadir/st31/reads/ont

    cat ~/data/20180827_2011_nina_fungus2/called/barcode04/*fastq.gz \
	~/data/20181108_0055_181107_nina_v2/called/barcode12/*fastq.gz \
	> $datadir/st90853/reads/ont/st90853.fastq.gz
    
    cat ~/data/20180827_2011_nina_fungus2/called/barcode05/*fastq.gz \
	~/data/20181108_0055_181107_nina_v2/called/barcode11/*fastq.gz \
	> $datadir/st31/reads/ont/st31.fastq.gz
    
fi

if [ $1 == longfq ] ; then
    for i in st31 st5317 st90853 ;
    do
	python3 ~/Code/utils/fastq_long.py \
		-i $datadir/$i/reads/ont/$i.fastq.gz \
		-o $datadir/$i/reads/ont/${i}_long.fastq.gz \
		-l 3000
    done
fi

if [ $1 == canu ] ; then
    ##for i in st31 st5317 st90853 ;
    for i in st90853 ;
    do
	mkdir -p $datadir/$i/asm/canu
	canu \
	    -p $i \
	    -d $datadir/$i/asm/canu \
	    genomeSize=39m \
	    -nanopore-raw $datadir/$i/reads/ont/${i}_long.fastq.gz
    done
fi


if [ $1 == flye ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/asm/flye
	flye \
	    --nano-raw $datadir/$i/reads/ont/${i}_long.fastq.gz \
	    -o $datadir/$i/asm/flye \
	    -g 39m \
	    -t 54
    done
fi


if [ $1 == trim_illumina ] ; then
    for i in st31 st5317 st90853 ;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	     -threads 36 -phred33 \
	     $datadir/$i/reads/illumina/$i*R1*.fastq.gz $datadir/$i/reads/illumina/$i*R2*.fastq.gz \
	     $datadir/$i/reads/illumina/${i}_fwd_paired.fq.gz $datadir/$i/reads/illumina/${i}_fwd_unpaired.fq.gz \
	     $datadir/$i/reads/illumina/${i}_rev_paired.fq.gz $datadir/$i/reads/illumina/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == nameflye ] ; then
    for i in st31 st5317 st90853 ;
    do
	cp $datadir/$i/asm/flye/assembly.fasta $datadir/$i/asm/flye/$i.ctgs.fasta
    done
fi
	
