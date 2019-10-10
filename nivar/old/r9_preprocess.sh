#!/bin/bash

##assemble on whack
prefix=nivar_r9
datadir=/kyber/Data/seqlab/sp_2019/$prefix

##called using guppy 3.2.1 in docker container

if [ $1 == fqs ] ; then
    cat $datadir/called/*fastq > $datadir/$prefix.fq
fi

    
if [ $1 == longfqs ] ; then
    python2 ~/Code/utils/fastq_long.py -i $datadir/$prefix.fq -o $datadir/${prefix}_3k.fq -l 3000
fi


if [ $1 == assemble ] ; then
    canu \
	-p nivar -d $datadir/canu \
	genomeSize=15m \
	-nanopore-raw $datadir/${prefix}_3k.fq
fi

if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed
    for i in cDNA gDNA ;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 72 -phred33 \
	     $datadir/${i}_illumina/CANI_${i}_R1.fastq.gz $datadir/${i}_illumina/CANI_${i}_R2.fastq.gz \
	     $datadir/trimmed/CANI_${i}_forward_paired.fq.gz $datadir/trimmed/CANI_${i}_forward_unpaired.fq.gz \
	     $datadir/trimmed/CANI_${i}_reverse_paired.fq.gz $datadir/trimmed/CANI_${i}_reverse_unpaired.fq.gz \
	     ILLUMINACLIP:idt_ud_indexes.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi


if [ $1 == pilon ] ; then
    mkdir -p $datadir/pilon
    cp $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz $datadir/pilon
    cp $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz $datadir/pilon

    bash ./pilon.sh $datadir/pilon $datadir/canu/nivar.contigs.fasta nivar
fi

if [ $1 == pilon_bwa ] ; then
    mkdir -p $datadir/pilon_bwa
    cp $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz $datadir/pilon_bwa
    cp $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz $datadir/pilon_bwa

    bash ./pilon_bwa.sh $datadir/pilon_bwa $datadir/canu/nivar.contigs.fasta nivar
fi

    
if [ $1 == racon ] ; then
    mkdir -p $datadir/racon

    ##racon is idiotic in its handling of filneames, so have to remove whitespace
    zcat $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz > $datadir/racon/CANI_gDNA_forward_paired.fq
    zcat $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz > $datadir/racon/CANI_gDNA_reverse_paired.fq

    sed -i -e 's/ //g' $datadir/racon/CANI_gDNA_forward_paired.fq
    sed -i -e 's/ //g' $datadir/racon/CANI_gDNA_reverse_paired.fq

    gzip $datadir/racon/*.fq
    bash ./racon.sh $datadir/racon $datadir/canu/nivar.contigs.fasta nivar
fi

    
if [ $1 == freebayes ] ; then
    mkdir -p $datadir/freebayes

    cp $datadir/canu/nivar.contigs.fasta $datadir/freebayes
    bowtie2-build $datadir/freebayes/nivar.contigs.fasta $datadir/freebayes/nivar.contigs
    bowtie2 -x $datadir/freebayes/nivar.contigs -1 $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz -2 $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz | samtools view -@ 24 -bS - | samtools sort -@ 24 -o $datadir/freebayes/nivar.sorted.bam
    samtools index $datadir/freebayes/nivar.sorted.bam
    
    freebayes -f $datadir/canu/nivar.contigs.fasta $datadir/freebayes/nivar.sorted.bam > $datadir/freebayes/nivar.vcf
fi
    
if [ $1 == freebayes_bwa ] ; then
    ##assumes you've already run the regular freebayes portion of stuff    
    mkdir -p $datadir/freebayes_bwa
    cp $datadir/trimmed/CANI_gDNA_forward_paired.fq.gz $datadir/freebayes_bwa
    cp $datadir/trimmed/CANI_gDNA_reverse_paired.fq.gz $datadir/freebayes_bwa

    bash ./freebayes_bwa.sh $datadir/freebayes_bwa $datadir/canu/nivar.contigs.fasta nivar
    
fi
    

if [ $1 == racon_bwa ] ; then
    ##add to racon analysis
    rm $datadir/racon/all.fq.gz
    bash ./racon_bwa.sh $datadir/racon $datadir/canu/nivar.contigs.fasta nivar_bwa
fi


if [ $1 == medusa ] ; then
    scadir=$datadir/medusa
    olddir=/kyber/Data/seqlab/sp_2019/fungus_asm
    mkdir -p $scadir

    cp -r ~/software/medusa/medusa_scripts ./

    java -jar ~/software/medusa/medusa.jar -f $olddir/References -i $datadir/canu/nivar.contigs.fasta -v -o $scadir/nivar_r10_scaffold.fasta
fi
