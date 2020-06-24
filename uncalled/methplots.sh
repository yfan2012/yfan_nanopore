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


if [ $1 == bisulfite ] ; then
    ##data from here unzipped:
    ##https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2308nnn/GSM2308632/suppl/GSM2308632_ENCFF279HCL_methylation_state_at_CpG_GRCh38.bed.gz
    ##https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2308nnn/GSM2308633/suppl/GSM2308633_ENCFF835NTC_methylation_state_at_CpG_GRCh38.bed.gz

    for i in $datadir/bisulfite/*.bed ;
    do
	prefix=`basename $i .bed | cut -d _ -f 2`
	bedtools intersect -wb -a $i -b $promos | \
	    awk '$14 =="promoter" {print $0}' > $datadir/bisulfite/promoters_$prefix.tsv
    done
fi
	    
