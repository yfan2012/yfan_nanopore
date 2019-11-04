#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis

for i in antibody infected mock Antibody_2dpi Antibody_3dpi Sindbis_2dpi Sindbis_3dpi ;
do
    if [ $1 == stringtie ] ; then
	mkdir -p $datadir/$i/stringtie
	stringtie -p 36 -L -G $datadir/rattus_norvegicus.gff -o $datadir/$i/stringtie/$i.rat.gtf $datadir/$i/align/$i.rat.splicealn.sorted.bam
    fi

    if [ $1 == transcriptaln ] ; then
	minimap2 -a -k14 -uf -t 36 $datadir/refs/rattus_norvegicus.rna.fa $datadir/$i/fqs/$i.fq \
	    | samtools view -@ 36 -b \
	    | samtools sort -@ 36 -o $datadir/$i/align/$i.rat.transcriptaln.sorted.bam -T $datadir/$i/align/reads.tmp -

	samtools index $datadir/$i/align/$i.rat.transcriptaln.sorted.bam

	samtools view -@ 36 -b -F 0x100 $datadir/$i/align/$i.rat.transcriptaln.sorted.bam \
	    | samtools sort -@ 36 -o $datadir/$i/align/$i.rat.transcriptaln.primary.sorted.bam
	samtools index $datadir/$i/align/$i.rat.transcriptaln.primary.sorted.bam
    fi
	
    
    if [ $1 == count_transcripts ] ; then
	mkdir -p $datadir/$i/transcript_counts
	python2 rat_expression.py \
	       -i $datadir/$i/align/$i.rat.transcriptaln.primary.sorted.bam \
	       -r $datadir/refs/rattus_norvegicus.rna.fa \
	       -o $datadir/$i/transcript_counts/$i.rat_transcript_counts.csv
    fi

done
       
