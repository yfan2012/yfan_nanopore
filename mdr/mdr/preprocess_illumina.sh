#!/bin/bash

prefix=181127_hiC_stool
datadir=/mithril/Data/Nanopore/projects/methbin/mdr
rawdir=/kyber/Data/NGS/Raw/${prefix}/FASTQ

if [ $1 == orgname ] ; then
    ##according to the csv in rawdir
    mkdir -p $datadir/illumina
    mkdir -p $datadir/illumina/raw
    cp $rawdir/HJFF7BCX2_2_AAGAGGCA~ACTGCATA_1.fastq.gz $datadir/illumina/raw/${prefix}_shotgun_1.fastq.gz
    cp $rawdir/HJFF7BCX2_2_AAGAGGCA~ACTGCATA_2.fastq.gz $datadir/illumina/raw/${prefix}_shotgun_2.fastq.gz
    cp $rawdir/HJFF7BCX2_2_TCGCCTTA~TAGATCGC_1.fastq.gz $datadir/illumina/raw/${prefix}_phase_1.fastq.gz
    cp $rawdir/HJFF7BCX2_2_TCGCCTTA~TAGATCGC_2.fastq.gz $datadir/illumina/raw/${prefix}_phase_2.fastq.gz
fi

if [ $1 == trim ] ; then
    ##trim
    untrimmed=$datadir/illumina/raw
    trimmed=$datadir/illumina/trimmed
    mkdir -p $trimmed
    cp ~/software/Trimmomatic-0.39/adapters/all_adapters.fa ./
    for i in shotgun phase ;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $untrimmed/${prefix}_${i}_1.fastq.gz $untrimmed/${prefix}_${i}_2.fastq.gz \
	     $trimmed/${prefix}_${i}_fwd_paired.fq.gz $trimmed/${prefix}_${i}_fwd_unpaired.fq.gz \
	     $trimmed/${prefix}_${i}_rev_paired.fq.gz $trimmed/${prefix}_${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == assemble ] ; then
    ##spades on shotgun reads. probs not worth assembling hic reads

    mkdir -p $datadir/illumina/metaspades
    cat $datadir/illumina/trimmed/*shotgun*unpaired* > $datadir/illumina/trimmed/${prefix}_shotgun_all_unpaired.fq.gz
    
    spades.py --meta \
	      --plasmid \
	      -t 36 \
	      -o $datadir/illumina/metaspades \
	      -1 $datadir/illumina/trimmed/${prefix}_shotgun_fwd_paired.fq.gz \
	      -2 $datadir/illumina/trimmed/${prefix}_shotgun_rev_paired.fq.gz \
	      -s $datadir/illumina/trimmed/${prefix}_shotgun_all_unpaired.fq.gz
    mv $datadir/illumina/metaspades/contigs.fasta $datadir/illumina/metaspades/${prefix}_shotgun.contigs.fasta
fi

spadesasm=$datadir/illumina/metaspades/${prefix}_shotgun.contigs.fasta
if [ $1 == illumina_amr ] ; then
    mkdir -p $datadir/illumina/amr
    for i in ecoh card ncbi resfinder plasmidfinder vfdb ecoli_vf megares argannot ;
    do
	abricate \
	    --threads 36 \
	    --db $i \
	    $spadesasm > $datadir/illumina/amr/$prefix.$i.tsv
    done
fi
   
if [ $1 == blast ] ; then
    mkdir -p $datadir/blast_contigs

    blastn \
	-num_threads 36 \
	-query $spadesasm \
	-db /atium/Data/ref/ncbi/nt \
	-outfmt 7 \
	-out $datadir/blast_contigs/$prefix.tsv
fi

