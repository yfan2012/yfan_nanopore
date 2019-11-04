#!/bin/bash

##Do some standard things for when you have illumina data
prefix=$1

rawdir=/uru/Data/NGS/Raw/$prefix
datadir=/uru/Data/dunlop/$prefix
refdir=/uru/Data/dunlop/reference

if [ $2 == trim ] ; then
    mkdir -p $datadir/trimmed
    cat ~/software/Trimmomatic-0.39/adapters/*.fa > ~/Code/yfan_nanopore/dunlop/all_adapters.fa
    for i in $rawdir/*_R1_*fastq.gz ;
    do
	samp=`basename $i .fastq.gz | cut -d _ -f 1`
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $i $rawdir/${samp}*_R2_*.fastq.gz \
	     $datadir/trimmed/${samp}_fwd_paired.fq.gz $datadir/trimmed/${samp}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${samp}_rev_paired.fq.gz $datadir/trimmed/${samp}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $2 == assemble ] ; then
    mkdir -p $datadir/spades
    for i in $datadir/trimmed/*_fwd_paired.fq.gz ;
    do
	samp=`basename $i _fwd_paired.fq.gz`
	mkdir -p $datadir/spades/$samp
	spades.py -1 $datadir/trimmed/${samp}_fwd_paired.fq.gz -2 $datadir/trimmed/${samp}_rev_paired.fq.gz -t 54 -m 300 -o $datadir/spades/$samp
	for file in $datadir/spades/$samp/* ;
	do
	    filename=`basename $file`
	    dir=`dirname $file`
	    mv $file $dir/${samp}_spades_$filename
	done
    done
fi

if [ $2 == parsnp ] ; then
    mkdir -p $datadir/parsnp
    mkdir -p $datadir/parsnp/genomes
    cp $datadir/spades/*/*_spades_contigs.fasta $datadir/parsnp/genomes/

    parsnp -p 36 -r ! -d $datadir/parsnp/genomes -o $datadir/parsnp -c
    harvesttools -i $datadir/parsnp/parsnp.ggr -V $datadir/parsnp/strain_snps.vcf
fi

if [ $2 == count_snps ] ; then
    Rscript ~/Code/utils/count_snps.R -i $datadir/parsnp/strain_snps.vcf -o $datadir/parsnp
fi

	
if [ $2 == align ] ; then
    mkdir -p $datadir/align
    ref=$refdir/bw_BW25113.fa

    bwa index $ref

    for i in $datadir/trimmed/*_fwd_paired.fq.gz ;
    do
	samp=`basename $i _fwd_paired.fq.gz`
	bwa mem -t 36 $ref $datadir/trimmed/${samp}_fwd_paired.fq.gz $datadir/trimmed/${samp}_rev_paired.fq.gz |\
	    samtools view -@ 36 -bS |\
	    samtools sort -@ 36 -o $datadir/align/$samp.sorted.bam
	samtools index $datadir/align/$samp.sorted.bam

    done
fi

if [ $2 == clean_bam ] ; then
    ##merge this into the align step if possible
    ##I think this crazy piping doesn't work because it needs the index of everyting besides the namesorted step 
    for i in $datadir/align/*.sorted.bam ;
    do
	samp=`basename $i .sorted.bam`
	samtools sort -@ 36 -n -o $datadir/align/$samp.namesorted.bam $i 
	
	samtools fixmate -@ 36 -r -m $datadir/align/$samp.namesorted.bam - |\
	    samtools sort -@ 36 -o $datadir/align/$samp.sorted.fxmt.bam
	samtools index $datadir/align/$samp.sorted.fxmt.bam
	
	samtools markdup -@ 36 -r $datadir/align/$samp.sorted.fxmt.bam - |\
	    samtools sort -@ 36 -o $datadir/align/$samp.sorted.mkdp.bam
	samtools index $datadir/align/$samp.sorted.mkdp.bam
    done
fi

if [ $2 == callvars ] ; then
    cd ~/software/freebayes/scripts
    mkdir -p $datadir/vars
   
    for i in $datadir/align/*.sorted.bam ;
    do
	samp=`basename $i .sorted.bam`
	
	##from freebayes docs, for naive reporting of snps and indels
	./freebayes-parallel <(./fasta_generate_regions.py $refdir/bw_BW25113.fa.fai 100000) 36 \
				   -f $refdir/bw_BW25113.fa \
				   --haplotype-length 0 \
				   --min-alternate-count 1 \
				   --min-alternate-fraction 0 \
				   --pooled-continuous \
				   $datadir/align/$samp.sorted.mkdp.bam > $datadir/vars/$samp.vcf
    done
fi

       
if [ $2 == varsupport ] ; then
    for i in $datadir/vars/*.vcf ;
    do
	samp=`basename $i .vcf`
	touch $datadir/vars/total_vars.csv
	python find_snps.py -v $i -g $refdir/bw_BW25113.gff3 -o $datadir/vars/$samp.csv >> $datadir/vars/total_vars.csv
    done
fi
