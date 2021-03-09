#!/bin/bash

##ragtag_fc identified as best asm for st31
##nina idenfieds contig22 and scaffold 19 as potentially one contig

datadir=/pym/Data/Nanopore/projects/prolificans/st31
asm=$datadir/ragtag/fc/st31.ragtag_fc.scaffolds.fasta
fq=$datadir/reads/ont/st31_long.fastq.gz

if [ $1 == align ] ; then
    mkdir -p $datadir/mancheck
    prefix=`basename $asm .fasta`
    
    minimap2 -t 36 -ax map-ont $asm $fq |
	samtools view -@ 36 -b |
	samtools sort -@ 36 -o $datadir/mancheck/$prefix.sorted.bam
    samtools index $datadir/mancheck/$prefix.sorted.bam
fi

if [ $1 == getpaf ] ; then
    ##paf file will probably be easier to work with for what I want
    prefix=`basename $asm .fasta`
    minimap2 -t 36 -x map-ont $asm $fq \
	     > $datadir/mancheck/$prefix.paf
fi

if [ $1 == filtpaf ] ; then
    prefix=`basename $asm .fasta`
    grep contig_22_RagTag $datadir/mancheck/$prefix.paf | \
	cut -d $'\t' -f 1-12 \
	    > $datadir/mancheck/$prefix.contig_22_RagTag.paf
    grep scaffold_19_RagTag $datadir/mancheck/$prefix.paf | \
	cut -d $'\t' -f 1-12 \
	    > $datadir/mancheck/$prefix.scaffold_19_RagTag.paf
fi
