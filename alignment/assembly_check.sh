#!/bin/bash

illreads=/atium/Data/NGS/Raw/171002_bucc
npreads=/atium/Data/Nanopore/cpowgs/170816_BUCC/170816_BUCC.fq
aligndir=/atium/Data/Nanopore/cpowgs/170816_BUCC/assembly_align

ref=$aligndir/170816_BUCC.pilon.fasta


##bowtie2-build $ref $aligndir/index/170816_BUCC.pilon

##bowtie2 -p 12 -x $aligndir/index/170816_BUCC.pilon -1 $illreads/bucc_S1_L001_R1_001.fastq.gz -2 $illreads/bucc_S1_L001_R2_001.fastq.gz | samtools view -bS | samtools sort -o $aligndir/170816_BUCC.ill.pilon.sorted.bam
##samtools index $aligndir/170816_BUCC.ill.pilon.sorted.bam

minimap2 -t 12 -a -x map-ont $ref $npreads | samtools view -b - | samtools sort -o $aligndir/170816_BUCC.np.pilon.sorted.bam
samtools index $aligndir/170816_BUCC.np.pilon.sorted.bam

