#!/bin/bash

datadir=~/work/r10_zymo

if [ $1 == data ] ; then
    mkdir -p $datadir/loman_asm
    ##wget http://nanopore.s3.climb.ac.uk/mockcommunity/v3/7cd60d3b-eafb-48d1-9aab-c8701232f2f8.ctg.cns.fa is r10 data supposedly
    ##wget http://nanopore.s3.climb.ac.uk/mockcommunity/v2/Zymo-GridION-EVEN-BB-SN.fq.wtdbg2.L5000-e3-K10000-NMX6000-S1-p23-r1.ctg.lay.gz.fa is r9 data supposedly
fi

if [ $1 == genomefix ] ; then
    for i in $datadir/Reference/zymo_refs/Genomes/*fasta ;
    do
	prefix=`basename $i | cut -d _ -f 1`
	sed -i -e "s/>/>$prefix/g" $i
    done
    cat $datadir/Reference/zymo_refs/Genomes/*fasta > $datadir/Reference/zymo_refs/all_genomes.fasta
fi


if [ $1 == mummer ] ; then
    mkdir -p $datadir/mummer

    nucmer -p $datadir/mummer/allgenomes_r10 $datadir/Reference/zymo_refs/all_genomes.fasta $datadir/loman_asm/7cd60d3b-eafb-48d1-9aab-c8701232f2f8.ctg.cns.fa
    nucmer -p $datadir/mummer/allgenomes_r9 $datadir/Reference/zymo_refs/all_genomes.fasta $datadir/loman_asm/Zymo-GridION-EVEN-BB-SN.fq.wtdbg2.L5000-e3-K10000-NMX6000-S1-p23-r1.ctg.lay.gz.fa
    ##mummerplot --png $datadir/mummer/allgenomes.layout $datadir/mummer/allgenomes.delta -R $datadir/Reference/zymo_refs/all_genomes.fasta -Q $datadir/loman_asm/7cd60d3b-eafb-48d1-9aab-c8701232f2f8.ctg.cns.fa
    ##mummerplot --png $datadir/mummer/allgenomes $datadir/mummer/allgenomes.delta

    ##dnadiff -p $datadir/mummer/allgenomes $datadir/Reference/zymo_refs/all_genomes.fasta $datadir/loman_asm/7cd60d3b-eafb-48d1-9aab-c8701232f2f8.ctg.cns.fa
    dnadiff -p $datadir/mummer/allgenomes_r10 -d $datadir/mummer/allgenomes_r10.delta
    dnadiff -p $datadir/mummer/allgenomes_r9 -d $datadir/mummer/allgenomes_r9.delta
fi

if [ $1 == diffs ] ; then
    mkdir -p $datadir/diffs

    python ~/Code/utils/motif_enrich.py -s $datadir/mummer/allgenomes_r10.snps -r $datadir/Reference/zymo_refs/all_genomes.fasta -m 6 -o $datadir/diffs/allgenomes_r10_6mer.csv
    python ~/Code/utils/motif_enrich.py -s $datadir/mummer/allgenomes_r10.snps -r $datadir/Reference/zymo_refs/all_genomes.fasta -m 5 -o $datadir/diffs/allgenomes_r10_5mer.csv
    python ~/Code/utils/motif_enrich.py -s $datadir/mummer/allgenomes_r10.snps -r $datadir/Reference/zymo_refs/all_genomes.fasta -m 4 -o $datadir/diffs/allgenomes_r10_4mer.csv

    python ~/Code/utils/motif_enrich.py -s $datadir/mummer/allgenomes_r9.snps -r $datadir/Reference/zymo_refs/all_genomes.fasta -m 6 -o $datadir/diffs/allgenomes_r9_6mer.csv
    python ~/Code/utils/motif_enrich.py -s $datadir/mummer/allgenomes_r9.snps -r $datadir/Reference/zymo_refs/all_genomes.fasta -m 5 -o $datadir/diffs/allgenomes_r9_5mer.csv
    python ~/Code/utils/motif_enrich.py -s $datadir/mummer/allgenomes_r9.snps -r $datadir/Reference/zymo_refs/all_genomes.fasta -m 4 -o $datadir/diffs/allgenomes_r9_4mer.csv
    
fi
