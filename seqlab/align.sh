#!/bin/bash

ml samtools
datadir=~/work/seqlab/fungus

if [ $1 == cDNA_canu ] ; then
    mkdir -p $datadir/align/canu
    
    cp $datadir/canu_pilon/candida_nivariensis_canu.pilon.6.fasta $datadir/align/canu
    sed -i -e 's/_pilon//g' $datadir/align/canu/candida_nivariensis_canu.pilon.6.fasta

    hisat2-build $datadir/align/canu/candida_nivariensis_canu.pilon.6.fasta $datadir/align/canu/candida_nivariensis_canu.pilon.6
    hisat2 -p 36 -q -x $datadir/align/canu/candida_nivariensis_canu.pilon.6.hisat -1 $datadir/trimmed/CANI_cDNA_forward_paired.fq.gz -2 $datadir/trimmed/CANI_cDNA_reverse_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/canu/CANI_canu_pilon_cDNA.sorted.bam
    samtools index $datadir/align/canu/CANI_canu_pilon_cDNA.sorted.bam
fi
   
if [ $1 == cDNA_spades ] ; then
    mkdir -p $datadir/align/spades
    cp $datadir/spades/contigs.fasta $datadir/align/spades/candida_nivariensis_spades.fasta
    hisat2-build $datadir/align/spades/candida_nivariensis_spades.fasta $datadir/align/spades/candida_nivariensis_spades
    hisat2 -p 36 -q -x $datadir/align/spades/candida_nivariensis_spades.hisat -1 $datadir/trimmed/CANI_cDNA_forward_paired.fq.gz -2 $datadir/trimmed/CANI_cDNA_reverse_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/spades/CANI_spades_cDNA.sorted.bam
    samtools index $datadir/align/spades/CANI_spades_cDNA.sorted.bam
fi

if [ $1 == cDNA_ref ] ; then
    mkdir -p $datadir/align/ref
    cp $datadir/References/candida_nivariensis.fa $datadir/align/ref
    hisat2-build $datadir/align/ref/candida_nivariensis.fa $datadir/align/ref/candida_nivariensis_ref
    hisat2 -p 36 -q -x $datadir/align/ref/candida_nivariensis_ref -1 $datadir/trimmed/CANI_cDNA_forward_paired.fq.gz -2 $datadir/trimmed/CANI_cDNA_reverse_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/ref/CANI_ref_cDNA.sorted.bam
    samtools index $datadir/align/ref/CANI_ref_cDNA.sorted.bam
fi

if [ $1 == dRNA_canu ] ; then
    mkdir -p $datadir/align/canu
    minimap2 -a -x splice -uf -k14 -t 36 $datadir/align/canu/candida_nivariensis_canu.pilon.6.fasta $datadir/fastqs/candida_nivariensis_dRNA.fastq | samtools view -b | samtools sort -o $datadir/align/canu/CANI_canu_pilon_dRNA.sorted.bam -T $datadir/align/canu/reads.tmp -
    samtools index $datadir/align/canu/CANI_canu_pilon_dRNA.sorted.bam
fi

if [ $1 == dRNA_spades ] ; then
    mkdir -p $datadir/align/spades
    minimap2 -a -x splice -uf -k14 -t 36 $datadir/align/spades/candida_nivariensis_spades.fasta $datadir/fastqs/candida_nivariensis_dRNA.fastq | samtools view -b | samtools sort -o $datadir/align/spades/CANI_spades_dRNA.sorted.bam -T $datadir/align/spades/reads.tmp -
    samtools index $datadir/align/spades/CANI_spades_dRNA.sorted.bam
fi

if [ $1 == dRNA_ref ] ; then
    mkdir -p $datadir/align/ref
    minimap2 -a -x splice -uf -k14 -t 36 $datadir/align/ref/candida_nivariensis.fa $datadir/fastqs/candida_nivariensis_dRNA.fastq | samtools view -b | samtools sort -o $datadir/align/ref/CANI_ref_dRNA.sorted.bam -T $datadir/align/ref/reads.tmp -
    samtools index $datadir/align/ref/CANI_ref_dRNA.sorted.bam
fi
