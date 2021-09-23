#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == trf ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/repeats
	asm=$datadir/$i/final/$i.final2.fasta
	
	~/software/trf $asm 2 7 7 80 10 50 600 -ngs > $datadir/$i/repeats/$i.trf.txt
    done
fi
	


if [ $1 == repeatmodel_db ] ; then
    mkdir -p ~/data/prolificans
    for i in st31 st5317 st90853 ;
    do
	mkdir -p ~/data/prolificans/$i
	cp $datadir/$i/final/$i.final2.fasta ~/data/prolificans/$i/
	~/software/RepeatModeler/BuildDatabase \
	    -name ~/data/prolificans/$i/$i.final2 \
	    ~/data/prolificans/$i/$i.final2.fasta
    done
fi

if [ $1 == repeatmodel ] ; then
    mkdir -p ~/data/prolificans
    for i in st31 st5317 st90853 ;
    do
	mkdir -p ~/data/prolificans/$i
	nohup ~/software/RepeatModeler/RepeatModeler \
	    -database ~/data/prolificans/$i/$i.final2 \
	    -pa 12 -LTRStruct >& ~/data/prolificans/$i/$i.out &
    done
    ##moved output dirs RM* manually to datadir
fi



if [ $1 == repeatmodel_test ] ; then
    mkdir -p ~/data/prolificans
    for i in st31 st5317 st90853 ;
    do
	mkdir -p ~/data/prolificans/$i
	nohup ~/software/RepeatModeler/RepeatModeler \
	    -database ~/data/prolificans/$i/$i.final2 \
	    -pa 12 >& ~/data/prolificans/$i/$i.out &
    done
fi



if [ $1 == repeatmasker ] ; then
    for i in st31 st5317 st90853 ;
    do
	~/software/RepeatMasker/RepeatMasker \
	    -pa 12 \
	    -s \
	    -species Fungi \
	    -gff \
	    $datadir/$i/final/$i.final2.fasta
    done
    
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/repeats/repeatmasker
	mv $datadir/$i/final/$i.final2.fasta.* $datadir/$i/repeats/repeatmasker/
    done
fi


if [ $1 == buildhmmlib ] ; then
    for i in st31 st5317 st90853 ;
    do
	sed -i -e '/^$/d' $datadir/$i/repeats/repeatmodeler/consensi.fa

	hmmbuild \
	    --dna \
	    --cpu 12 \
	    $datadir/$i/repeats/repeatmodeler/consensi.hmm \
	    $datadir/$i/repeats/repeatmodeler/families.stk
    done
fi

if [ $1 == classifyreps ] ; then
    #for i in st31 st5317 st90853 ;
    for i in st5317 st90853 ;
    do
	~/software/RepeatModeler/RepeatClassifier \
	    -consensi $datadir/$i/repeats/repeatmodeler/consensi.fa \
	    -stockholm $datadir/$i/repeats/repeatmodeler/families.stk
    done
fi

    
if [ $1 == repeatmasker_custom ] ; then
    for i in st31 st5317 st90853 ;
    do
	mkdir -p $datadir/$i/repeats/repeatmasker_custom
	cp $datadir/$i/final/$i.final2.fasta $datadir/$i/repeats/repeatmasker_custom
	~/software/RepeatMasker/RepeatMasker \
	    -pa 12 \
	    -s \
	    -lib $datadir/$i/repeats/repeatmodeler/consensi.fa.classified \
	    -gff \
	    $datadir/$i/repeats/repeatmasker_custom/$i.final2.fasta
    done
fi
