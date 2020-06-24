#!/bin/bash

datadir=/archive/skovaka/nanopore_runs
ref=/archive/skovaka/refs/human/hg38_noalt.fa

if [ $1 == combinefq ] ; then
	cat $datadir/20191101_GM12878_control/combined.fastq \
		$datadir/20191212_GM12878_control/combined.fq \
		$datadir/20191220_GM12878_control/combined.fq \
		$datadir/inv_all/rebc_flank20k/combined.fq > ~/data/fastqs/ontwgs.fq
	cat $datadir/20191101_GM12878_control/combined.fastq \
		$datadir/20191212_GM12878_control/combined.fq \
		$datadir/20191220_GM12878_control/combined.fq > ~/data/fastqs/control.fq
fi	

if [ $1 == index ] ; then
	cp $datadir/inv_all/rebc_flank20k/bc/sequencing_summary.txt ~/data/ontwgs_seqsum.txt
	tail -n +2 $datadir/20191101_GM12878_control/knope_20191101_131853_FAL59666_minion_sequencing_run_2019_GM12878_control_sequencing_summary.txt >> ~/data/ontwgs_seqsum.txt
	tail -n +2 $datadir/20191212_GM12878_control/bc/sequencing_summary.txt >> ~/data/ontwgs_seqsum.txt
	tail -n +2 $datadir/20191220_GM12878_control/bc/sequencing_summary.txt >> ~/data/ontwgs_seqsum.txt
	~/software/nanopolish/nanopolish index \
		-d $datadir/inv_all/rebc_flank20k/fast5s \
		-d $datadir/20191101_GM12878_control \
		-d $datadir/20191212_GM12878_control \
		-d $datadir/20191220_GM12878_control \
		-s ~/data/ontwgs_seqsum.txt \
		~/data/fastqs/ontwgs.fq



	##grab data we sequenced just incase
	cp $datadir/20191101_GM12878_control/knope_20191101_131853_FAL59666_minion_sequencing_run_2019_GM12878_control_sequencing_summary.txt ~/data/control_seqsum.txt
	tail -n +2 $datadir/20191212_GM12878_control/bc/sequencing_summary.txt >> ~/data/control_seqsum.txt
	tail -n +2 $datadir/20191220_GM12878_control/bc/sequencing_summary.txt >> ~/data/control_seqsum.txt
	~/software/nanopolish/nanopolish index \
		-d $datadir/20191101_GM12878_control \
		-d $datadir/20191212_GM12878_control \
		-d $datadir/20191220_GM12878_control \
		-s ~/data/control_seqsum.txt \
		~/data/fastqs/control.fq
fi	

if [ $1 == align ] ; then
	for i in ontwgs control ; 
	do
		~/software/minimap2/minimap2 -ax map-ont -t 16 /archive/skovaka/refs/human/hg38_noalt.mmi ~/data/fastqs/$i.fq |\
			samtools view -@ 15 -bS |\
			samtools sort -@ 15 -o ~/data/align/$i.sorted.bam
	done	
fi

if [ $1 == polish ] ; then
	mkdir -p ~/data/methcalls
	regions=/archive/skovaka/refs/genes/inv_all/genes_flank20k.bed
	for i in ontwgs control ;
	do	
		while read r; 
		do 
			chr=`echo $r | tr -s ' ' | cut -d ' ' -f 1`
			first=`echo $r | tr -s ' ' | cut -d ' ' -f 2`
			last=`echo $r | tr -s ' ' | cut -d ' ' -f 3`
			name=`echo $r | tr -s ' ' | cut -d ' ' -f 4`
			window=${chr}:${first}-${last}
			~/software/nanopolish/nanopolish call-methylation -t 16 \
				-r ~/data/fastqs/$i.fq \
				-b ~/data/align/$i.sorted.bam \
				-g $ref -w $window >> ~/data/methcalls/${i}_methcalls.tsv 
		done < $regions
	done
fi	


if [ $1 == methfreq ] ; then
        for i in ontwgs control ; 
	do	
		sed -i '/num_motifs/d' ~/data/methcalls/${i}_methcalls.tsv
        	header=`head -n 1 ~/data/methcalls/methcalls.tsv`
        	sed -i "1i $header" ~/data/methcalls/${i}_methcalls.tsv
        	~/software/nanopolish/scripts/calculate_methylation_frequency.py ~/data/methcalls/${i}_methcalls.tsv > ~/data/methcalls/${i}_methcalls_freq.tsv
	done
fi

