#!/bin/bash

datadir=/uru/Data/Nanopore/projects/uncalled
promos=/uru/Data/Nanopore/projects/uncalled/ref/hg38_promoters.gff

if [ $1 == picksites ] ; then
    for i in $datadir/methcalls/*_freq.tsv ;
    do
	prefix=`basename $i .tsv | cut -d _ -f 1`
	tail --lines=+2 $i | \
	    bedtools intersect -wb -a stdin -b $promos | \
	    awk '$11 =="promoter" {print $0}' > $datadir/methcalls/promoters_$prefix.tsv
    done
fi

if [ $1 == phasemeth ] ; then
    outdir=~/Dropbox/timplab_data/uncalled/methview
    vcf=$datadir/ref/hg38.hybrid.vcf
    bed=$datadir/ref/genes.bed

    echo INDEXING+++++++++++++++++++++++++++++++++++++++++++
    bedtools intersect -a $vcf -b $bed > $datadir/ref/hg38_targeted.hybrid.vcf
    ##whatshap giving me a weird error. remove the line for now
    sed -i '/GTTTTTGGGTTGCTGG/d' $datadir/ref/hg38_targeted.hybrid.vcf
    
    bgzip $datadir/ref/hg38_targeted.hybrid.vcf
    tabix -p vcf $datadir/ref/hg38_targeted.hybrid.vcf.gz

    echo HAPLOTAG+++++++++++++++++++++++++++++++++++++++++++
    for i in $outdir/*colored.sorted.bam ;
    do
	prefix=`basename $i .bam`
	samtools index $i
	whatshap haplotag --ignore-read-groups \
		 -o $datadir/align/$prefix.haplotag.bam \
		 $datadir/ref/hg38_targeted.hybrid.vcf.gz  \
		 $i
	samtools index $datadir/align/$prefix.haplotag.bam
    done
fi


if [ $1 == splitbam ] ; then
    for i in $datadir/align/*.haplotag.bam ;
    do
	prefix=`basename $i .haplotag.bam`
	python v2_split_BAM_by_tag.py \
	       -b $datadir/align/$prefix.haplotag.bam \
	       -o $datadir/align \
	       -n $prefix
	   
	   samtools index $datadir/align/*hap1.bam
	   samtools index $datadir/align/*hap2.bam
    done
fi

