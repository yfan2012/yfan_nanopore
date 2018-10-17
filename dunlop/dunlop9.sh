#!/bin/bash

datadir=~/work/181004_dunlop_9ecoli
srcdir=~/Code/utils/marcc

mkdir -p $datadir/batch_logs

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_dunlop $srcdir/untar.scr $datadir/181004_dunlop_9ecoli.tar.gz $datadir
fi

   
if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-365 --job-name=call_dunlop --output=$datadir/call_logs/dunlop_ecoli3.%A_%a.out $srcdir/bc_call.scr $datadir
fi

if [ $1 == rmraw ] ; then
    for i in $datadir/raw/* ;
    do
	(rm -r $i) &
    done
fi

if [ $1 == fastq ] ; then
    mkdir -p $datadir/fastqs
    for i in 1 2 3 4 5 6 7 8 9 ;
    do
	cat $datadir/called/*/workspace/pass/barcode0${i}/*fastq > $datadir/fastqs/ecoli${i}.fastq
    done
fi
    
if [ $1 == assemble ] ; then
    ##everything uses canu 1.7 by default now
    mkdir -p $datadir/assemblies
    for i in 1 2 3 4 5 6 7 8 9 ;
    do
	mkdir $datadir/assemblies/ecoli${i}
	bash assemble_bacteria.sh $datadir/fastqs/ecoli${i}.fastq $datadir/assemblies/ecoli${i}
    done
fi

    
if [ $1 == illumina_marcc ] ; then
    mkdir -p ~/work/181004_dunlop_9ecoli/illumina
    scp -r smaug:/dilithium/Data/NGS/Raw/181013_dunlop9_GEOB/ecoli* ~/work/181004_dunlop_9ecoli/illumina/
fi

if [ $1 == reorg ] ; then
    for i in 1 2 3 4 5 6 7 8 9 ;
    do
	mkdir -p $datadir/ecoli${i}
	mv $datadir/assemblies/ecoli${i} $datadir/ecoli${i}/assembly
	mv $datadir/fastqs/ecoli${i}.fastq $datadir/ecoli${i}/
	mkdir -p $datadir/ecoli${i}/pilon
	mv $datadir/illumina/ecoli9-${i}* $datadir/ecoli${i}/pilon/
    done
fi
	

if [ $1 == pilon ] ; then
    for i in 1 2 3 4 5 6 7 8 9 ;
    do
	mkdir -p $datadir/ecoli${i}/pilon
	sbatch --output=$datadir/batch_logs/ecoli${i}_pilon.out --job-name=ecoli${i} pilon9.scr ecoli${i}
    done
fi




if [ $1 == copyassembly ] ; then
    scp -r $datadir/ecoli* smaug:/dilithium/Data/Nanopore/projects/dunlop_ecoli/
fi


if [ $1 == parsnp ] ; then
    parsdir=~/Dropbox/yfan/dunlop/parsnp/ecoli9
    assembledir=~/Dropbox/yfan/dunlop/assemblies
    
    mkdir -p $parsdir

    parsnp -r $assembledir/bw25311.fasta -d $assembledir/ecoli9 -p 12 -o $parsdir -c
    harvesttools -i $parsdir/parsnp.ggr -V $parsdir/raw.vcf
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
