#!/bin/bash

##picking up from years ago

datadir=/atium/Data/projects/jose_burkle
asm=$datadir/followup_plasmid/Jose3.fasta

if [ $1 == clearcarriage ] ; then
    sed -i -e 's/\r//g' $asm
fi


    
if [ $1 == plasmid ] ; then
    ##follow up on sus segment that might not be plasmid
    ##seq is from jose email

    abricate \
	--threads 36 \
	--db plasmidfinder \
	$asm > $datadir/followup_plasmid/Jose3.tsv
fi


if [ $1 == mummer_self ] ; then
   ##no results from plasmid finder at all - see if p3 is a repeated contig

   mkdir -p $datadir/followup_plasmid/mummer_self
   
   nucmer \
       -p $datadir/followup_plasmid/mummer_self/Jose3 \
       $asm \
       $asm

   mummerplot --postscript \
	      -p $datadir/followup_plasmid/mummer_self/Jose3 \
	      $datadir/followup_plasmid/mummer_self/Jose3.delta

   dnadiff \
       -p $datadir/followup_plasmid/mummer_self/Jose3 \
       $asm \
       $asm

fi

if [ $1 == showcoords ] ; then
    show-coords -l -T $datadir/followup_plasmid/mummer_self/Jose3.delta \
		> $datadir/followup_plasmid/mummer_self/Jose3.sc.tsv
fi

		
if [ $1 == mum_self ] ; then
    mkdir -p $datadir/followup_plasmid/mum
    mummer \
	-mum \
	-l 500 \
	-b \
	-n \
	-qthreads 12 \
	$asm \
	$asm > $datadir/followup_plasmid/mum/selfaligns.out
fi

    
if [ $1 == align ] ; then
    mkdir -p $datadir/followup_plasmid/align

    minimap2 -t 36 -ax map-ont $asm $datadir/np_fastqs/jose3.fq |
	samtools view -@ 36 -b |
	samtools sort -@ 36 -o $datadir/followup_plasmid/align/Jose3.sorted.bam
    samtools index $datadir/followup_plasmid/align/Jose3.sorted.bam
fi

	     
	
    
