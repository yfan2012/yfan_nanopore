#!/bin/bash

datadir=/kyber/Data/Nanopore/projects/tellurite
dbox=~/Dropbox/Timplab_Data/tellurite
abriout=$dbox/abricate
blastout=$dbox/blast_ter
repout=$dbox/report
prokout=$dbox/prokka

if [ $1 == get_ter_genes ] ; then
    ##get all ter c and d seqs using trick from geo
    esearch -db protein -query 'txid573[orgn] AND terC[prot]' | efetch -format fasta_cds_na > $datadir/ref/cds.terC.fa
    esearch -db protein -query 'txid573[orgn] AND terD[prot]' | efetch -format fasta_cds_na > $datadir/ref/cds.terD.fa
fi


if [ $1 == plasfind ] ; then
    for qual in raw pilon ;
    do
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    abricate --db plasmidfinder --quiet $i > $abriout/$qual/$prefix.plasmidfinder.tsv
	done
    done
fi

if [ $1 == mkblastdb ] ; then
    for qual in raw pilon ;
    do
	mkdir -p $datadir/blast_db/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    makeblastdb -in $i -out $datadir/blast_db/$qual/$prefix.db -dbtype nucl
	done
    done
fi

if [ $1 == search_ter ] ; then
    for qual in raw pilon ;
    do
	mkdir -p $blastout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
    	    prefix=`basename $i .fasta`
	    blastn -query $datadir/ref/cds.terC.fa -db $datadir/blast_db/$qual/$prefix.db -outfmt 7 -out $blastout/$qual/$prefix.terC.tsv
	    sed -i -e 's/_pilon//g' $blastout/$qual/$prefix.terC.tsv
	    blastn -query $datadir/ref/cds.terD.fa -db $datadir/blast_db/$qual/$prefix.db -outfmt 7 -out $blastout/$qual/$prefix.terD.tsv
	    sed -i -e 's/_pilon//g' $blastout/$qual/$prefix.terD.tsv
	done
    done
fi

if [ $1 == report ] ; then
    for qual in raw pilon ;
    do
	mkdir -p $repout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    echo $prefix
	    sed -i -e 's/_pilon//g' $abriout/$qual/$prefix.plasmidfinder.tsv
	    python ~/Code/yfan_nanopore/tellurite/tellurite_plasfind.py -b $blastout/$qual/$prefix.terD.tsv -p $abriout/$qual/$prefix.plasmidfinder.tsv -g terD -n $prefix -o $repout/$qual/$prefix.terD.report.csv
	    python ~/Code/yfan_nanopore/tellurite/tellurite_plasfind.py -b $blastout/$qual/$prefix.terC.tsv -p $abriout/$qual/$prefix.plasmidfinder.tsv -g terC -n $prefix -o $repout/$qual/$prefix.terC.report.csv
	done
    done
fi

if [ $1 == prokka ] ; then
    ##for qual in raw pilon ;
    for qual in pilon ;
    do
	mkdir -p $prokout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    prokka --outdir $prokout/$qual --genus Klebsiella --usegenus --prefix $prefix --cpus 12 --force $i &> $prokout/$qual/$prefix.log.txt
	done
    done
fi


