#!/bin/bash

datadir=~/work/180714_dunlop_3ecoli
srcdir=~/Code/utils/marcc

mkdir -p $datadir/batch_logs

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_dunlop $srcdir/untar.scr $datadir/180714_dunlop_3ecoli.tar.gz $datadir
fi

   
if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-693 --job-name=call_dunlop --output=$datadir/call_logs/dunlop_ecoli3.%A_%a.out $srcdir/bc_call.scr $datadir
fi


if [ $1 == fastq ] ; then
    mkdir -p $datadir/fastqs
    cat $datadir/called/*/workspace/pass/barcode01/*fastq > $datadir/fastqs/ecoli1.fastq
    cat $datadir/called/*/workspace/pass/barcode02/*fastq > $datadir/fastqs/ecoli2.fastq
    cat $datadir/called/*/workspace/pass/barcode03/*fastq > $datadir/fastqs/ecoli3.fastq
fi
    
if [ $1 == downsamp ] ; then
    head -800000 $datadir/fastqs/ecoli1.fastq > $datadir/fastqs/ecoli1_sub200k.fastq
    head -800000 $datadir/fastqs/ecoli2.fastq > $datadir/fastqs/ecoli2_sub200k.fastq
    head -800000 $datadir/fastqs/ecoli3.fastq > $datadir/fastqs/ecoli3_sub200k.fastq
fi

if [ $1 == assemble ] ; then
    mkdir $datadir/ecoli1_assembly
    bash assemble_bacteria.sh $datadir/fastqs/ecoli1.fastq $datadir/ecoli1_assembly
    mkdir $datadir/ecoli2_assembly
    bash assemble_bacteria.sh $datadir/fastqs/ecoli2.fastq $datadir/ecoli2_assembly
    mkdir $datadir/ecoli3_assembly
    bash assemble_bacteria.sh $datadir/fastqs/ecoli3.fastq $datadir/ecoli3_assembly
fi

if [ $1 == assemble17 ] ; then
    mkdir $datadir/ecoli1_assembly17
    bash assemble_bacteria.sh $datadir/fastqs/ecoli1_sub200k.fastq $datadir/ecoli1_assembly17
    mkdir $datadir/ecoli2_assembly17
    bash assemble_bacteria.sh $datadir/fastqs/ecoli2_sub200k.fastq $datadir/ecoli2_assembly17
    mkdir $datadir/ecoli3_assembly17
    bash assemble_bacteria.sh $datadir/fastqs/ecoli3_sub200k.fastq $datadir/ecoli3_assembly17
fi

    
if [ $1 == illumina_dilith ] ; then
    ##copy illumina data to dilithium
    mkdir -p /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli
    for i in `find ~/BaseSpace/Projects/180722_dunlop_3ecoli/ -name *fastq.gz` ;
    do
	cp $i /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/
    done
fi

if [ $1 == illumina_marcc ] ; then
    mkdir -p ~/work/180714_dunlop_3ecoli/illumina
    scp -r smaug:/dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/* ~/work/180714_dunlop_3ecoli/illumina/
fi

if [ $1 == pilon ] ; then
    sbatch --output=$datadir/batch_logs/ecoli1_pilon.out --job-name=ecoli1 pilon.scr ecoli1
    sbatch --output=$datadir/batch_logs/ecoli2_pilon.out --job-name=ecoli2 pilon.scr ecoli2
    sbatch --output=$datadir/batch_logs/ecoli3_pilon.out --job-name=ecoli3 pilon.scr ecoli3
fi


if [ $1 == copyassembly ] ; then
    scp -r $datadir/ecoli* smaug:/dilithium/Data/Nanopore/projects/dunlop_ecoli/
fi


if [ $1 == parsnp ] ; then
    parsdir=~/Dropbox/yfan/dunlop/parsnp
    assembledir=~/Dropbox/yfan/dunlop/assemblies
    
    mkdir -p $parsdir/raw_parsnp
    mkdir -p $parsdir/canu17_parsnp
    mkdir -p $parsdir/pilon_parsnp
    mkdir -p $parsdir/pilon17_parsnp

    parsnp -r $assembledir/bw25311.fasta -d $assembledir/raw -p 12 -o $parsdir/raw_parsnp -c
    harvesttools -i $parsdir/raw_parsnp/parsnp.ggr -V $parsdir/raw_parsnp/raw.vcf
    parsnp -r $assembledir/bw25311.fasta -d $assembledir/canu17 -p 12 -o $parsdir/canu17_parsnp -c
    harvesttools -i $parsdir/canu17_parsnp/parsnp.ggr -V $parsdir/canu17_parsnp/canu17.vcf
    parsnp -r $assembledir/bw25311.fasta -d $assembledir/pilon -p 12 -o $parsdir/pilon_parsnp -c
    harvesttools -i $parsdir/pilon_parsnp/parsnp.ggr -V $parsdir/pilon_parsnp/pilon.vcf
    parsnp -r $assembledir/bw25311.fasta -d $assembledir/pilon17 -p 12 -o $parsdir/pilon17_parsnp -c
    harvesttools -i $parsdir/pilon17_parsnp/parsnp.ggr -V $parsdir/pilon17_parsnp/pilon17.vcf
fi

if [ $1 == spades ] ; then
    ml python/2.7
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	fastqdir=~/work/180714_dunlop_3ecoli/illumina
	sampdir=$datadir/${i}_spades
	mkdir $sampdir
	mkdir $fastqdir/trimmed
	
	java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $fastqdir/$i*R1* $fastqdir/$i*R2* $fastqdir/trimmed/${i}_forward_paired.fq.gz $fastqdir/trimmed/${i}_forward_unpaired.fq.gz $fastqdir/trimmed/${i}_reverse_paired.fq.gz $fastqdir/trimmed/${i}_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	
	spades.py -1 $fastqdir/trimmed/${i}_forward_paired.fq.gz -2 $fastqdir/trimmed/${i}_reverse_paired.fq.gz -t 36 -m 300 -o $datadir/${i}_spades
    done
fi



if [ $1 == reorg ] ; then
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	mv -p $datadir/$i
	mv $datadir/${i}_* $datadir/$i
	mv $datadir/$i/*assembly $datadir/$i/assembly
	mv $datadir/$i/*assembly17 $datadir/$i/assembly17
	mv $datadir/$i/*pilon $datadir/$i/pilon
	mv $datadir/$i/*pilon17 $datadir/$i/pilon17
    done
fi
	
if [ $1 == align_pilon17diff ] ; then
    ##look at persistent CCWGG snp in ecoli 3 by alignment
    ml samtools
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	fastqdir=~/work/180714_dunlop_3ecoli/illumina
	bamdir=$datadir/$i/align
	pilonill=$datadir/$i/align/pilonill
	
	mkdir -p $bamdir
	mkdir -p $pilonill
	mkdir -p $pilonill/btidx
	
	cp $datadir/$i/pilon17/$i.pilon17.10.fasta $pilonill/btidx/
	bowtie2-build -q $pilonill/btidx/$i.pilon17.10.fasta $pilonill/btidx/ecoli3.pilon17.10
	bowtie2 -p 24 -x $pilonill/btidx/ecoli3.pilon17.10 -1 $datadir/illumina/$i*R1_001.fastq.gz -2 $datadir/illumina/$i*R2_001.fastq.gz | samtools view -bS - | samtools sort -o $pilonill/$i.sorted.bam
	samtools index $pilonill/$i.sorted.bam
    done
fi

	
if [ $1 == cppilonill ] ; then
   for i in ecoli1 ecoli2 ecoli3 ;
   do
       pilonill=$datadir/$i/align/pilonill
       scp -r $pilonill/$i.sorted.bam* smaug:~/Dropbox/yfan/dunlop/align/pilonill/
   done
fi
       
if [ $1 == align_pilondiff ] ; then
    ##look at persistent CCWGG snp in ecoli 3 by alignment
    ml samtools
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	fastqdir=~/work/180714_dunlop_3ecoli/illumina
	bamdir=$datadir/$i/align
	pilonill=$datadir/$i/align/pilonill
	
	mkdir -p $bamdir
	mkdir -p $pilonill
	mkdir -p $pilonill/btidx
	
	cp $datadir/$i/pilon/$i.pilon.10.fasta $pilonill/btidx/
	bowtie2-build -q $pilonill/btidx/$i.pilon.10.fasta $pilonill/btidx/ecoli3.pilon.10
	bowtie2 -p 24 -x $pilonill/btidx/ecoli3.pilon.10 -1 $datadir/illumina/$i*R1_001.fastq.gz -2 $datadir/illumina/$i*R2_001.fastq.gz | samtools view -bS - | samtools sort -o $pilonill/$i.sorted.bam
	samtools index $pilonill/$i.sorted.bam
    done
fi
