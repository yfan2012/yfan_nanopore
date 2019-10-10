#!/bin/bash

##datadir=~/Dropbox/timplab_data/birdpair
datadir=/kyber/Data/Nanopore/projects/birdpair

if [ $1 == cleanfiles ] ; then
    sed -i -e 's/_pilon//g' $datadir/bird_pilon/bird_canu.pilon.10*
    sed -i -e 's/>/>pilon_/g' $datadir/bird_pilon/bird_canu.pilon.10*
    sed -i -e 's/_pilon//g' $datadir/patient_pilon/patient_canu.pilon.10*
    sed -i -e 's/>/>pilon_/g' $datadir/patient_pilon/patient_canu.pilon.10*
fi


if [ $1 == mummerfilt ] ; then
    ##do all the mummer things
    mkdir -p $datadir/mummer_filter
    mkdir -p ~/tmp

    cp $datadir/bird_pilon/bird_canu.pilon.10.fasta ~/tmp/
    cp $datadir/patient_pilon/patient_canu.pilon.10.fasta ~/tmp/
    cp $datadir/ref/CRNE_$i
    
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


if [ $1 == cleanmwrsa852 ] ; then
    ##parsnp bitches if there are dashes in the seq name soooo
    fafile=$datadir/ref/CRNE_mwrsa852.fa
    sed -i -e 's/-//g' $fafile
fi

for i in mwrsa852 c45 grubii ; do
    if [ $1 == mummer ] ; then
	##do all the mummer things
	mkdir -p $datadir/mummer_$i
	mkdir -p ~/tmp

	cp $datadir/bird_pilon/bird_canu.pilon.10.fasta ~/tmp/
	cp $datadir/patient_pilon/patient_canu.pilon.10.fasta ~/tmp/
	cp $datadir/ref/CRNE_$i.fa ~/tmp/

	echo ======================== running nucmer =========================
	nucmer -p ~/tmp/bird_$i ~/tmp/CRNE_$i.fa ~/tmp/bird_canu.pilon.10.fasta
	nucmer -p ~/tmp/patient_$i ~/tmp/CRNE_$i.fa ~/tmp/patient_canu.pilon.10.fasta

	echo ======================== running mummerplot 2 ===================
	mummerplot --fat --png -p ~/tmp/bird_$i ~/tmp/bird_$i.delta
	mummerplot --fat --png -p ~/tmp/patient_$i ~/tmp/patient_$i.delta

	echo ======================== running dnadiff =========================
	dnadiff -p ~/tmp/bird_$i ~/tmp/CRNE_$i.fa ~/tmp/bird_canu.pilon.10.fasta
	dnadiff -p ~/tmp/patient_$i ~/tmp/CRNE_$i.fa ~/tmp/patient_canu.pilon.10.fasta

	
	echo ======================== cleaning up ==============================
	mkdir -p $datadir/mummer/bird_$i
	mkdir -p $datadir/mummer/patient_$i
	cp ~/tmp/bird_$i* $datadir/mummer/bird_$i/
	cp ~/tmp/patient_$i* $datadir/mummer/patient_$i/
	rm -r ~/tmp/*
    fi

    if [ $1 == parsnp ] ; then
	mkdir -p $datadir/parsnp_$i
	mkdir -p $datadir/parsnp_$i/genomes
	
	cp $datadir/*spades/*spades.fasta $datadir/parsnp_$i/genomes
	cp $datadir/ref/CRNE_$i.fa $datadir/parsnp_$i/genomes
	
	parsnp -p 36 -r $datadir/parsnp_$i/genomes/CRNE_$i.fa -d $datadir/parsnp_$i/genomes -o $datadir/parsnp_$i -c
	harvesttools -i $datadir/parsnp_$i/parsnp.ggr -V $datadir/parsnp_$i/parsnp.vcf
    fi
    
    if [ $1 == count_snps ] ; then
	out=~/Dropbox/timplab_data/birdpair/snp_annot_$i.csv
	gff=$datadir/ref/CRNE_$i.gff
	vcf=$datadir/parsnp_$i/parsnp.vcf
	
	python ~/Code/yfan_nanopore/birdpair/find_snps.py -v $vcf -g $gff -o $out
    fi
    
    if [ $1 == count_unique ] ; then
	outfile=~/Dropbox/timplab_data/birdpair/snp_annot_$i.csv
	##great little awk thing that counts number of lines that are unique by first two columns
	echo Number of unique coding snps in $i : 
	awk -F "," '!seen[$1,$2]++' $outfile | wc -l
    fi

    if [ $1 == filt_snps ] ; then
	awk '$7 == "PASS" {print}' $datadir/parsnp_$i/parsnp.vcf > $datadir/parsnp_$i/parsnp_filt.vcf
    fi
    
    
    if [ $1 == count_filt_snps ] ; then
	out=~/Dropbox/timplab_data/birdpair/snp_annot_filt_$i.csv
	gff=$datadir/ref/CRNE_$i.gff
	vcf=$datadir/parsnp_$i/parsnp_filt.vcf
	
	python ~/Code/yfan_nanopore/birdpair/find_snps.py -v $vcf -g $gff -o $out
    fi
    
    if [ $1 == count_filt_unique ] ; then
	out=~/Dropbox/timplab_data/birdpair/snp_annot_filt_$i.csv
	##great little awk thing that counts number of lines that are unique by first two columns
	echo Number of unique coding filtered snps in $i :
	awk -F "," '!seen[$1,$2]++' $out | wc -l
    fi
done
