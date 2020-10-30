#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=/uru/Data/Nanopore/projects/nivar/paperfigs
fq=/uru/Data/Nanopore/projects/nivar/r9/r9_3kb.fq
dbxdir=~/Dropbox/yfan/nivar/paperfigs

if [ $1 == runqc ] ; then
    Rscript ~/Code/utils/qc/run_summary.R \
	    -i $rawdir/r9/run1_called/sequencing_summary.txt \
	    -o $dbxdir/run1_run_summary.pdf \
	    -p Run1
    Rscript ~/Code/utils/qc/run_summary.R \
	    -i $rawdir/r9/run2_called/sequencing_summary.txt \
	    -o $dbxdir/run2_run_summary.pdf \
	    -p Run2
    Rscript ~/Code/utils/qc/run_summary.R \
	    -i $rawdir/r10/called/sequencing_summary.txt \
	    -o $dbxdir/r10_summary.pdf \
	    -p r10
fi

if [ $1 == assemble ] ; then
    mkdir $datadir/assembly
    
    canu \
	-p nivar \
	-d $datadir/assembly \
	genomeSize=11.6m \
	-nanopore-raw $fq
fi
asm=$datadir/assembly/nivar.contigs.fasta

if [ $1 == medaka ] ; then
    mkdir -p $datadir/medaka
    
    medaka_consensus \
	-i $fq \
	-d $asm \
	-o $datadir/medaka \
	-t 36 \
	-m r941_min_high_g303
fi


if [ $1 == polish ] ; then
    mkdir -p $datadir/freebayes
    cp $rawdir/illumina/gDNA_trimmed/*_paired.fq.gz $datadir/freebayes
    
    bash ./freebayes_bwa.sh $datadir/freebayes $asm nivar
fi
asmcorr=$datadir/freebayes/nivar_fb3_bwa.fasta


ref=$rawdir/reference/candida_nivariensis.fa
gla=$rawdir/reference/medusa_fungi/candida_glabrata.fa


if [ $1 == glabrata_mum ] ; then
    mkdir -p ~/tmp/mummer
    mkdir -p ~/tmp/mummer/glabrata_nivar_fb3_bwa
    
    cp $gla ~/tmp/mummer

    cp $asmcorr ~/tmp/mummer/glabrata_nivar_fb3_bwa
    nucmer -p ~/tmp/mummer/glabrata_nivar_fb3_bwa/glabrata_nivar_fb3_bwa ~/tmp/mummer/glabrata_nivar_fb3_bwa/nivar_fb3_bwa.fasta ~/tmp/mummer/candida_glabrata.fa 
    mummerplot --filter --fat --png -p ~/tmp/mummer/glabrata_nivar_fb3_bwa/glabrata_nivar_fb3_bwa ~/tmp/mummer/glabrata_nivar_fb3_bwa/glabrata_nivar_fb3_bwa.delta
    dnadiff -p ~/tmp/mummer/glabrata_nivar_fb3_bwa/glabrata_nivar_fb3_bwa ~/tmp/mummer/glabrata_nivar_fb3_bwa/nivar_fb3_bwa.fasta ~/tmp/mummer/candida_glabrata.fa 

    cp -r ~/tmp/mummer $datadir/
fi

if [ $1 == align ] ; then
    mkdir $datadir/align

    minimap2 -t 36 -ax map-ont $asmcorr $fq |
	samtools view -@ 36 -b |
	samtools sort -@ 36 -o $datadir/align/nivar_fb3_bwa.sorted.bam
    samtools index $datadir/align/nivar_fb3_bwa.sorted.bam
fi


if [ $1 == mito_mum ] ; then
    mkdir -p ~/tmp/mummer/mito

    cp $datadir/assembly_final/nivar_fb3_bwa_mito.fasta ~/tmp/mummer/mito
    cp $rawdir/reference/mitos/candida_nivariensis_mito.fa ~/tmp/mummer/mito/
    
    nucmer -p ~/tmp/mummer/mito/nivar_mito ~/tmp/mummer/mito/nivar_fb3_bwa_mito.fasta ~/tmp/mummer/mito/candida_nivariensis_mito.fa 
    mummerplot --filter --fat --png -p ~/tmp/mummer/mito/nivar_mito ~/tmp/mummer/mito/nivar_mito.delta
    dnadiff -p ~/tmp/mummer/mito/nivar_mito ~/tmp/mummer/mito/nivar_fb3_bwa_mito.fasta ~/tmp/mummer/mito/candida_nivariensis_mito.fa
    cp -r ~/tmp/mummer/mito $datadir/mummer/
fi

fin=$datadir/assembly_final/nivar.final.fasta
if [ $1 == clean_names ] ; then
    sed -i -e 's/\s.*$//' $fin
    sed -i -e 's/00//g' $fin
fi
    
if [ $1 == medusa ] ; then
    mkdir -p $datadir/medusa

    cp -r ~/software/medusa/medusa_scripts ./

    mkdir -p $datadir/medusa/nivariensis
    java -jar ~/software/medusa/medusa.jar \
	 -f $rawdir/reference/medusa_fungi \
	 -i $fin \
	 -v $ref \
	 -o $datadir/medusa/nivariensis/nivar.final.scaffold.fasta

    mkdir -p $datadir/medusa/glabrata
    java -jar ~/software/medusa/medusa.jar \
	 -f $rawdir/reference/medusa_fungi \
	 -i $fin \
	 -v $gla \
	 -o $datadir/medusa/glabrata/nivar.final.scaffold.fasta
	
fi

if [ $1 == ragtag ] ; then
    mkdir -p $datadir/ragtag
    mkdir -p $datadir/ragtag/nivariensis
    ragtag.py scaffold \
	      -w \
	      -u \
	      -o $datadir/ragtag/nivariensis \
	      $ref \
	      $fin
    mkdir -p $datadir/ragtag/glabrata
    ragtag.py scaffold \
	      -w \
	      -u \
	      -o $datadir/ragtag/glabrata \
	      $gla \
	      $fin
	
fi

sca=$datadir/medusa/nivar.scaffold.fasta
rag=$datadir/ragtag/ragtag.scaffolds.fasta

if [ $1 == mummer ] ; then
    mkdir -p $datadir/mummer
    mkdir -p ~/tmp/mummer
    mkdir -p ~/tmp/mummer/scaffold_medusa
    mkdir -p ~/tmp/mummer/scaffold_ragtag
    mkdir -p ~/tmp/mummer/nivar_fb3_bwa
    
    cp $ref ~/tmp/mummer

    cp $sca ~/tmp/mummer/scaffold_medusa
    nucmer -p ~/tmp/mummer/scaffold_medusa/nivar.scaffold ~/tmp/mummer/scaffold_medusa/nivar.scaffold.fasta ~/tmp/mummer/candida_nivariensis.fa 
    mummerplot --filter --fat --png -p ~/tmp/mummer/scaffold_medusa/nivar.scaffold ~/tmp/mummer/scaffold_medusa/nivar.scaffold.delta
    dnadiff -p ~/tmp/mummer/scaffold_medusa/nivar.scaffold ~/tmp/mummer/scaffold_medusa/nivar.scaffold.fasta ~/tmp/mummer/candida_nivariensis.fa 

    cp $rag ~/tmp/mummer/scaffold_ragtag
    nucmer -p ~/tmp/mummer/scaffold_ragtag/ragtag.scaffolds ~/tmp/mummer/scaffold_ragtag/ragtag.scaffolds.fasta ~/tmp/mummer/candida_nivariensis.fa 
    mummerplot --filter --fat --png -p ~/tmp/mummer/scaffold_ragtag/ragtag.scaffolds ~/tmp/mummer/scaffold_ragtag/ragtag.scaffolds.delta
    dnadiff -p ~/tmp/mummer/scaffold_ragtag/ragtag.scaffolds ~/tmp/mummer/scaffold_ragtag/ragtag.scaffolds.fasta ~/tmp/mummer/candida_nivariensis.fa 

    cp $asmcorr ~/tmp/mummer/nivar_fb3_bwa
    nucmer -p ~/tmp/mummer/nivar_fb3_bwa/nivar_fb3_bwa ~/tmp/mummer/nivar_fb3_bwa/nivar_fb3_bwa.fasta ~/tmp/mummer/candida_nivariensis.fa 
    mummerplot --filter --fat --png -p ~/tmp/mummer/nivar_fb3_bwa/nivar_fb3_bwa ~/tmp/mummer/nivar_fb3_bwa/nivar_fb3_bwa.delta
    dnadiff -p ~/tmp/mummer/nivar_fb3_bwa/nivar_fb3_bwa ~/tmp/mummer/nivar_fb3_bwa/nivar_fb3_bwa.fasta ~/tmp/mummer/candida_nivariensis.fa 

    cp -r ~/tmp/mummer $datadir/
fi

