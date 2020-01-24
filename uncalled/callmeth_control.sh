#!/bin/bash

datadir=/archive/skovaka/nanopore_runs
ref=/archive/skovaka/refs/human/hg38_noalt.fa

if [ $1 == gatherfq ] ; then
	##get 20191101_GM12878_control
	cp $datadir/20191101_GM12878_control/combined.fastq ~/data/fastqs/20191101_GM12878_control.fq
	##get 20191212_GM12878_control
	cp $datadir/20191212_GM12878_control/combined.fq ~/data/fastqs/20191212_GM12878_control.fq
	##get 20191220_GM12878_control
	cp $datadir/20191220_GM12878_control/combined.fq ~/data/fastqs/20191220_GM12878_control.fq
	##get NA12878
	cp $datadir/NA12878/all.fastq ~/data/fastqs/NA12878.fq
fi

if [ $1 == combinefq ] ; then
	cat $datadir/inv_all/rebc/bc/*fastq > ~/data/fastqs/control.fq

fi	

if [ $1 == index ] ; then
	#cat $datadir/NA12878/rel_6_sequencing_summary.txt $datadir/20191220_GM12878_control/seqsum*.txt $datadir/20191212_GM12878_control/seqsum*.txt $datadir/20191101_GM12878_control/*sequencing_summary.txt > ~/data/control_seqsum.txt
	~/software/nanopolish/nanopolish index -d $datadir/inv_all/fast5s -s $datadir/inv_all/rebc/bc/sequencing_summary.txt ~/data/fastqs/control.fq
fi	

if [ $1 == align ] ; then
	samtools sort -@ 15 -o ~/data/align/control.sorted.bam -T tmp $datadir/inv_all/rebc/mm2_genes.bam
fi

if [ $1 == polish ] ; then
	mkdir -p ~/data/methcalls
	regions=/archive/skovaka/refs/genes/inv_all/genes_flank20k.bed
	rm ~/data/methcalls/control_methcalls.tsv
	while read r; 
	do 
		chr=`echo $r | tr -s ' ' | cut -d ' ' -f 1`
		first=`echo $r | tr -s ' ' | cut -d ' ' -f 2`
		last=`echo $r | tr -s ' ' | cut -d ' ' -f 3`
		name=`echo $r | tr -s ' ' | cut -d ' ' -f 4`
		window=${chr}:${first}-${last}
		~/software/nanopolish/nanopolish call-methylation -t 16 -r ~/data/fastqs/control.fq -b ~/data/align/control.sorted.bam -g $ref -w $window >> ~/data/methcalls/control_methcalls.tsv 
	done < $regions
fi	


if [ $1 == methfreq ] ; then
        sed -i '/num_motifs/d' ~/data/methcalls/control_methcalls.tsv
        header=`head -n 1 ~/data/methcalls/methcalls.tsv`
        sed -i "1i $header" ~/data/methcalls/control_methcalls.tsv
        ~/software/nanopolish/scripts/calculate_methylation_frequency.py ~/data/methcalls/control_methcalls.tsv > ~/data/methcalls/control_methcalls_freq.tsv
fi

