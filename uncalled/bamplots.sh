#!/bin/bash

##taken from how Isac and Tim do methylation plots in IGV
##uses Isac's python utils

datadir=/uru/Data/Nanopore/projects/uncalled


if [ $1 == mtsv2bed ] ; then
    ##nanopolish tsv to bedfile

    for i in methcalls control_methcalls ontwgs_methcalls ;
    do
	./mtsv2bedGraph.py -i $datadir/methcalls/$i.tsv > $datadir/methcalls/$i.bed
    done
fi

if [ $1 == index ] ; then
    ##index bedfile

    for i in methcalls control_methcalls ontwgs_methcalls ;
    do
	bedtools sort -i $datadir/methcalls/$i.bed > $datadir/methcalls/$i.sorted.bed
	bgzip $datadir/methcalls/$i.sorted.bed
	tabix -p bed $datadir/methcalls/$i.sorted.bed.gz
    done
fi


if [ $1 == bamcolor ] ; then
    ./convert_bam_for_methylation.py \
	-t 12 \
	-f $datadir/ref/hg38_noalt.fa \
	-b $datadir/align/20191220_GM12878_inv_allRU.sorted.bam \
	-c $datadir/methcalls/methcalls.sorted.bed.gz \
	-r $datadir/genes_flank20k.bed \
	-o $datadir/methcalls/methcalls_colored.bam
    ./convert_bam_for_methylation.py \
	-t 12 \
	-f $datadir/ref/hg38_noalt.fa \
	-b $datadir/align/control.sorted.bam \
	-c $datadir/methcalls/control_methcalls.sorted.bed.gz \
	-r $datadir/genes_flank20k.bed \
	-o $datadir/methcalls/control_methcalls_colored.bam
    ./convert_bam_for_methylation.py \
	-t 12 \
	-f $datadir/ref/hg38_noalt.fa \
	-b $datadir/align/ontwgs.sorted.bam \
	-c $datadir/methcalls/ontwgs_methcalls.sorted.bed.gz \
	-r $datadir/genes_flank20k.bed \
	-o $datadir/methcalls/ontwgs_methcalls_colored.bam

fi

	
if [ $1 == sort ] ; then
    samtools sort -@ 36 -o $datadir/methcalls/methcalls_colored.sorted.bam  $datadir/methcalls/methcalls_colored.bam
    samtools index $datadir/methcalls/methcalls_colored.sorted.bam
    samtools sort -@ 36 -o $datadir/methcalls/control_methcalls_colored.sorted.bam $datadir/methcalls/control_methcalls_colored.bam
    samtools index $datadir/methcalls/control_methcalls_colored.sorted.bam
    samtools sort -@ 36 -o $datadir/methcalls/ontwgs_methcalls_colored.sorted.bam $datadir/methcalls/ontwgs_methcalls_colored.bam
    samtools index $datadir/methcalls/ontwgs_methcalls_colored.sorted.bam
fi

