#!/bin/bash

##datadir=~/Dropbox/timplab_data/birdpair
datadir=/kyber/Data/Nanopore/projects/birdpair
asmdir=$datadir/assemblies/canu

if [ $1 == cleanfiles ] ; then
    sed -i -e 's/_pilon//g' $asmdir/bird_canu.pilon.10*
    sed -i -e 's/>/>pilon_/g' $asmdir/bird_canu.pilon.10*
    sed -i -e 's/_pilon//g' $asmdir/patient_canu.pilon.9*
    sed -i -e 's/>/>pilon_/g' $asmdir/patient_canu.pilon.9*
fi

if [ $1 == mummer ] ; then
    ##do all the mummer things
    mkdir -p $datadir/mummer
    mkdir -p ~/tmp

    cp $asmdir/bird_canu.pilon.10.fasta ~/tmp/
    cp $asmdir/patient_canu.pilon.9.fasta ~/tmp/
    
    nucmer -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    mummerplot --fat --png -p ~/tmp/birdpair.layout ~/tmp/birdpair.delta -R ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta
    mummerplot --fat --png -p ~/tmp/birdpair ~/tmp/birdpair.delta

    dnadiff -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    cp ~/tmp/birdpair* $datadir/mummer/
fi

if [ $1 == mummerfilt ] ; then
    ##do all the mummer things
    mkdir -p $datadir/mummer_filter
    mkdir -p ~/tmp

    cp $asmdir/bird_canu.pilon.10.fasta ~/tmp/
    cp $asmdir/patient_canu.pilon.9.fasta ~/tmp/
    
    nucmer -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    mummerplot --filter --fat --png -p ~/tmp/birdpair.layout ~/tmp/birdpair.delta -R ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta
    mummerplot --filter --fat --png -p ~/tmp/birdpair ~/tmp/birdpair.delta

    dnadiff -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    cp ~/tmp/birdpair* $datadir/mummer_filter/
fi


if [ $1 == trim ] ; then
    illdir=/kyber/Data/NGS/Raw/190214_nina_birdpair
    mkdir -p $datadir/trimmed
    cat ~/software/Trimmomatic-0.39/adapters/*fa > ~/Code/yfan_nanopore/birdpair/all.fa

    for i in bird patient ; do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 \
	     $illdir/$i*_R1_*.fastq.gz \
	     $illdir/$i*_R1_*.fastq.gz \
	     $datadir/trimmed/${i}_fwd_paired.fq.gz \
	     $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${i}_rev_paired.fq.gz \
	     $datadir/trimmed/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == spades ] ; then
    for i in bird patient ; do
	mkdir -p $datadir/${i}_spades
	spades.py -1 $datadir/trimmed/${i}_fwd_paired.fq.gz -2 $datadir/trimmed/${i}_rev_paired.fq.gz -t 54 -m 300 -o $datadir/${i}_spades
	mv $datadir/${i}_spades/contigs.fasta $datadir/${i}_spades/$i.spades.fasta
    done
fi


if [ $1 == parsnp ] ; then
    mkdir -p $datadir/parsnp
    mkdir -p $datadir/parsnp/genomes

    cp $datadir/*spades/*spades.fasta $datadir/parsnp/genomes
    cp $datadir/ref/CRNE_grubii.fa $datadir/parsnp/genomes

    parsnp -p 24 -r $datadir/parsnp/genomes/CRNE_grubii.fa -d $datadir/parsnp/genomes -o $datadir/parsnp -c
    harvesttools -i $datadir/parsnp/parsnp.ggr -V $datadir/parsnp/parsnp.vcf

fi

	     
	 
    
    
