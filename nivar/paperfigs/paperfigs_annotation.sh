#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=$rawdir/paperfigs
fq=$rawdir/r9/r9_3kb.fq
drna=$rawdir/dRNA/dRNA.fq
rnafwd=$rawdir/illumina/cDNA_trimmed/nivar_cDNA_fwd_paired.fq
rnarev=$rawdir/illumina/cDNA_trimmed/nivar_cDNA_rev_paired.fq
dbxdir=~/Dropbox/yfan/nivar/paperfigs

rawname=$datadir/assembly_final/nivar.final.raw_tignames.fasta
gen=$datadir/assembly_final/nivar.final.fasta
ref=$rawdir/reference/candida_nivariensis.fa

gla=$rawdir/reference/medusa_fungi/candida_glabrata.fa
cer=$rawdir/reference/saccharomyces_cerevisiae.fa
alb=$rawdir/reference/candida_albicans.fa
glagff=$rawdir/reference/candida_glabrata.gff
cergff=$rawdir/reference/saccharomyces_cerevisiae.gff
albgff=$rawdir/reference/candida_albicans.gff

if [ $1 == fix_tignames ] ; then
    mv $gen $rawname
    awk -F ' ' '{print $1}' > $gen
fi

if [ $1 == liftoff ] ; then
    mkdir -p $datadir/annotation/liftoff
    liftoff \
	-p 36 \
	-o $datadir/annotation/liftoff/nivar_cer_lifted.gff \
	-u $datadir/annotation/liftoff/nivar_cer_unmapped.txt \
	-dir $datadir/annotation/liftoff/cer_intermediate_files \
	-g $cergff \
	$gen \
	$cer

    liftoff \
	-p 36 \
	-o $datadir/annotation/liftoff/nivar_alb_lifted.gff \
	-u $datadir/annotation/liftoff/nivar_alb_unmapped.txt \
	-dir $datadir/annotation/liftoff/alb_intermediate_files \
	-g $albgff \
	$gen \
	$alb

    liftoff \
	-p 36 \
	-o $datadir/annotation/liftoff/nivar_gla_lifted.gff \
	-u $datadir/annotation/liftoff/nivar_gla_unmapped.txt \
	-dir $datadir/annotation/liftoff/gla_intermediate_files \
	-g $glagff \
	$gen \
	$gla

fi

if [ $1 == align ] ; then
    mkdir -p $datadir/annotation/align

    minimap2 -t 36 \
	-ax splice -uf -k14 $gen $drna | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/annotation/align/nivar_dRNA.sorted.bam
    samtools index $datadir/annotation/align/nivar_dRNA.sorted.bam

    echo hisat####################################
    hisat2-build $gen $datadir/assembly_final/nivar.final
    hisat2 -p 36 \
	-x $datadir/assembly_final/nivar.final \
	-1 $rnafwd \
	-2 $rnarev | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/annotation/align/nivar_cDNA.sorted.bam
    samtools index $datadir/annotation/align/nivar_cDNA.sorted.bam
fi

if [ $1 == stringtie_drna ] ; then
    mkdir -p $datadir/annotation/stringtie

    stringtie \
	$datadir/annotation/align/nivar_dRNA.sorted.bam \
	-o $datadir/annotation/stringtie/denovo_drna.gff \
	-L \
	-p 36

fi

if [ $1 == stringtie_rnaseq ] ; then
    mkdir -p $datadir/annotation/stringtie

    stringtie \
	$datadir/annotation/align/nivar_cDNA.sorted.bam \
	-o $datadir/annotation/stringtie/denovo_rnaseq.gff \
	-p 36
fi

if [ $1 == braker ] ; then
    mkdir -p $datadir/annotation/braker

    export GENEMARK_PATH=~/software/gmes_linux_64
    ~/software/BRAKER/scripts/braker.pl \
	--cores=36 \
	--gff3 \
	--genome=$gen \
	--species=nivar \
	--bam=$datadir/annotation/align/nivar_cDNA.sorted.bam
fi
    
if [ $1 == gffcompare ] ; then
    mkdir -p $datadir/annotation/compare

    liftcer=$datadir/annotation/liftoff/nivar_cer_lifted.gff
    liftalb=$datadir/annotation/liftoff/nivar_alb_lifted.gff
    liftgla=$datadir/annotation/liftoff/nivar_gla_lifted.gff
    braker=$datadir/annotation/braker/braker.gff3
    drna=$datadir/annotation/stringtie/denovo_drna.gff
    rnaseq=$datadir/annotation/stringtie/denovo_rnaseq.gff
    final=$datadir/annotation_final/nivar.final.gff

    ##awk '$3=="u" {print $4}' liftgla_liftcer.nivar_cer_lifted.gff.tmap | grep -f - ../liftoff/nivar_cer_lifted.gff

    cp $braker $final

    gffcompare -r $final -o $datadir/annotation/compare/liftgla_final $liftgla
    awk '$3 == "u" {print $4}' $datadir/annotation/liftoff/liftgla_final.nivar_gla_lifted.gff.tmap \
	| grep -f - $liftgla >> $final

    gffcompare -r $final -o $datadir/annotation/compare/liftcer_final $liftcer
    awk '$3 == "u" {print $4}' $datadir/annotation/liftoff/liftcer_final.nivar_cer_lifted.gff.tmap \
	| grep -f - $liftcer >> $final

    gffcompare -r $final -o $datadir/annotation/compare/liftalb_final $liftalb
    awk '$3 == "u" {print $4}' $datadir/annotation/liftoff/liftalb_final.nivar_alb_lifted.gff.tmap \
	| grep -f - $liftalb >> $final
    
    gffcompare -r $final -o $datadir/annotation/compare/drna_final $drna
    awk '$3 == "u" {print $4}' $datadir/annotation/stringtie/drna_final.denovo_drna.gff.tmap \
	| grep -f - $drna >> $final
    
    ##gffcompare -r $final -o $datadir/annotation/compare/braker_final $braker
    ##awk '$3 == "u" {print $4}' $datadir/annotation/braker/braker_final.braker.gff3.tmap \
	##| grep -f - $braker >> $final


    ##using above awk lines to count how many features, genes, exons are contributed each time
fi


if [ $1 == transcriptome_oldattempt_not_used ] ; then
    mkdir -p $datadir/annotation/combined
    Rscript ~/Code/yfan_nanopore/nivar/paperfigs/annot_compare.R

    	##-r $datadir/annotation/combined/lifted_all.gff \
    gffcompare \
	-r $datadir/annotation/combined/lifted_all_gla.gff \
	-o $datadir/annotation/combined/lifted_vs_data \
	$datadir/annotation/combined/data_all.gff

    mkdir -p $datadir/annotation_final
    Rscript ~/Code/yfan_nanopore/nivar/paperfigs/annot_combine.R
fi
