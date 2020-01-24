#!/bin/bash

datadir=/archive/skovaka/nanopore_runs/20191220_GM12878_inv_allRU/
ref=/archive/skovaka/refs/human/hg38_noalt.fa

if [ $1 == gatherfq ] ; then
	cat $datadir/bc/*fastq > ~/data/fastqs/20191220_GM12878_inv_allRU.fq
fi

if [ $1 == index ] ; then
	~/software/nanopolish/nanopolish index -d $datadir -s $datadir/bc/sequencing_summary.txt ~/data/fastqs/20191220_GM12878_inv_allRU.fq
fi

if [ $1 == align ] ; then
	ref=/archive/skovaka/refs/human/hg38_noalt.fa
	mkdir -p ~/data/align
	~/software/minimap2/minimap2 -ax map-ont -t 8 $ref ~/data/fastqs/20191220_GM12878_inv_allRU.fq | samtools view -@ 4 -bS | samtools sort -@ 4 -o ~/data/align/20191220_GM12878_inv_allRU.sorted.bam -T tmp 
fi

if [ $1 == polish ] ; then
	mkdir -p ~/data/methcalls
	regions=/archive/skovaka/refs/genes/inv_all/genes_flank20k.bed
	rm ~/data/methcalls/methcalls.tsv
	while read r; 
	do 
		chr=`echo $r | tr -s ' ' | cut -d ' ' -f 1`
	        first=`echo $r | tr -s ' ' | cut -d ' ' -f 2`
	        last=`echo $r | tr -s ' ' | cut -d ' ' -f 3`
	        name=`echo $r | tr -s ' ' | cut -d ' ' -f 4`
	        window=${chr}:${first}-${last}
	        ~/software/nanopolish/nanopolish call-methylation -t 16 -r ~/data/fastqs/20191220_GM12878_inv_allRU.fq -b ~/data/align/20191220_GM12878_inv_allRU.sorted.bam -g $ref -w $window >> ~/data/methcalls/methcalls.tsv 
        done < $regions
fi      

if [ $1 == methfreq ] ; then
	sed -i '/num_motifs/d' ~/data/methcalls/methcalls.tsv
	header=`head -n 1 ~/data/methcalls/control_methcalls.tsv`
	sed -i "1i $header" ~/data/methcalls/methcalls.tsv
	~/software/nanopolish/scripts/calculate_methylation_frequency.py ~/data/methcalls/methcalls.tsv > ~/data/methcalls/methcalls_freq.tsv
fi	
