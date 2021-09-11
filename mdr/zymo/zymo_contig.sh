#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/zymo/raw
ssddir=~/data/mdr/zymo
prefix=20190809_zymo_control
datadir=/mithril/Data/Nanopore/projects/methbin/zymo


##ref from https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa

if [ $1 == assemble ] ; then
    mkdir -p $datadir/flye

    flye \
	--nano-raw $datadir/fastq/$prefix/$prefix.fq.gz \
	-o $datadir/flye/$prefix \
	-t 36 \
	-g 100m \
	--meta
    mv $datadir/flye/$prefix/assembly.fasta $datadir/flye/$prefix/$prefix.fasta
    ##was rerun with v2.9
fi

fq=$datadir/fastq/$prefix/$prefix.fq.gz

if [ $1 == racon_align ] ; then
    mkdir -p $datadir/racon

    minimap2 -t 36 -ax map-ont $datadir/flye/$prefix/$prefix.fasta $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $datadir/racon/$prefix.sorted.bam
    samtools index $datadir/racon/$prefix.sorted.bam
fi

if [ $1 == racon ] ; then		   
    samtools view -@ 36 $datadir/racon/$prefix.sorted.bam > $datadir/racon/$prefix.sorted.sam
    racon -m 8 -x -6 -g -8 -w 500 -t 54\
	  $fq \
	  $datadir/racon/$prefix.sorted.sam \
	  $datadir/flye/$prefix/$prefix.fasta > $datadir/racon/$prefix.racon.contigs.fasta
fi


if [ $1 == medaka ] ; then
    mkdir -p $datadir/medaka

    medaka_consensus \
	-i $fq \
	-d $datadir/racon/$prefix.racon.contigs.fasta \
	-o $datadir/medaka \
	-t 54 \
	-m r941_min_high_g360
fi


if [ $1 == untar ] ; then
    mkdir -p $ssddir/raw/$prefix
    tar -xzf $rawdir/$prefix.tar.gz -C $ssddir/raw/$prefix
fi


if [ $1 == megalodon_polished ] ; then
    mkdir -p $ssddir/megalodon

    megalodon \
        $ssddir/raw/$prefix/$prefix/fast5 \
        --overwrite \
        --guppy-server-path "/usr/bin/guppy_basecall_server" \
        --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --reference $datadir/medaka/consensus.fasta \
        --outputs per_read_mods \
        --output-directory $ssddir/megalodon/${prefix}_polished \
        --write-mods-text \
        --devices "cuda:0" \
        --processes 36
fi

if [ $1 == move_polished ] ; then
    mkdir -p $datadir/megalodon/${prefix}_polished

    mv $ssddir/megalodon/${prefix}_polished/* $datadir/megalodon/${prefix}_polished/
fi

if [ $1 == megaidx_polished ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt.idx
fi

if [ $1 == barcode_polished ] ; then
    mkdir -p $datadir/barcode/${prefix}_polished
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt.idx \
           -r $datadir/medaka/consensus.fasta \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
           -o $datadir/barcode/${prefix}_polished/${prefix}_contigs_barcodes15.txt \
           -t 12 ;} &> $datadir/barcode/${prefix}_polished/${prefix}_contigs_barcodes15_time.txt
fi

if [ $1 == barcode_polished_50 ] ; then
    mkdir -p $datadir/barcode/${prefix}_polished
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/${prefix}_polished/per_read_modified_base_calls.txt.idx \
           -r $datadir/medaka/consensus.fasta \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
           -o $datadir/barcode/${prefix}_polished/${prefix}_contigs_barcodes50.txt \
           -t 36 ;} &> $datadir/barcode/${prefix}_polished/${prefix}_contigs_barcodes50_time.txt
fi


if [ $1 == align_polished ] ; then
    minimap2 -t 36 -x map-ont $datadir/medaka/consensus.fasta $datadir/fastq/$prefix/$prefix.fq.gz \
	     > $datadir/align/${prefix}_polished.paf
    
    minimap2 -t 36 -ax map-ont $datadir/medaka/consensus.fasta $datadir/fastq/$prefix/$prefix.fq.gz | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/align/${prefix}_polished.sorted.bam
    samtools index $datadir/align/${prefix}_polished.sorted.bam
fi

if [ $1 == motifcounts_polished ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -r $datadir/medaka/consensus.fasta \
	   -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
	   -a $datadir/align/${prefix}_polished.paf \
	   -m $datadir/barcode/${prefix}_polished/${prefix}_barcodes15.txt \
	   -o $datadir/barcode/${prefix}_polished/${prefix}_barcodes15_motifcounts.txt \
	   -q 40
fi

if [ $1 == motifcounts_polished_50 ] ; then
    python ~/Code/yfan_meth/utils/megalodon_barcode_filter.py \
	   -r $datadir/medaka/consensus.fasta \
	   -b ~/Code/yfan_nanopore/mdr/rebase/barcodes50.txt \
	   -a $datadir/align/${prefix}_polished.paf \
	   -m $datadir/barcode/${prefix}_polished/${prefix}_barcodes50.txt \
	   -o $datadir/barcode/${prefix}_polished/${prefix}_barcodes50_motifcounts.txt \
	   -q 40
fi


if [ $1 == makeblastdb ] ; then
    name=`echo $ref | cut -d . -f 1`
    makeblastdb \
	-in $ref \
	-out $name \
	-dbtype nucl
fi

if [ $1 == blast_polished ] ; then
    name=`echo $ref | cut -d . -f 1`
    blastn \
	-num_threads 36 \
	-query $datadir/medaka/consensus.fasta \
	-db $name \
	-outfmt 7 \
	-out $datadir/barcode/${prefix}_polished/${prefix}_polished_blast.tsv
fi

if [ $1 == mummer_check_asm ] ; then
    ref=/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa
    asm=$datadir/medaka/consensus.fasta
    
    mkdir -p $datadir/mummer

    nucmer \
	-p $datadir/mummer/polished_ref \
	$asm \
	$ref
    
    mummerplot --filter --fat --png \
	       -p $datadir/mummer/polished_ref \
	       $datadir/mummer/polished_ref.delta

    dnadiff \
	-p $datadir/mummer/polished_ref \
	$asm \
	$ref
    
fi

if [ $1 == mummer_listeria ] ; then
    ##from asm and blast, it looks like listeria and faecalis have a 35kb region like exactly in common
    refdir=/uru/Data/Nanopore/projects/read_class/zymo/ref/D6322.refseq/Genomes
    listeria=$refdir/Listeria_monocytogenes_complete_genome.fasta
    faecalis=$refdir/Enterococcus_faecalis_complete_genome.fasta

    nucmer \
	-p $datadir/mummer/listeria_faecalis \
	$listeria \
	$faecalis
    
    mummerplot --filter --fat --png \
	       -p $datadir/mummer/listeria_faecalis \
	       $datadir/mummer/listeria_faecalis.delta
    
    dnadiff \
	-p $datadir/mummer/listeria_faecalis \
	$listeria \
	$faecalis
fi

if [ $1 == id_contigs ] ; then
    Rscript contig_id.R
fi

if [ $1 == align_polished ] ; then
    mkdir -p $datadir/align_polished

    minimap2 -t 36 -ax map-ont $datadir/medaka/consensus_bacteria.fasta $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $datadir/align_polished/$prefix.sorted.bam
    samtools index $datadir/align_polished/$prefix.sorted.bam

    minimap2 -t 36 -x map-ont $datadir/medaka/consensus_bacteria.fasta $fq \
	     > $datadir/align_polished/$prefix.paf
fi

    

################################################################
asm=$datadir/flye/$prefix/$prefix.fasta
if [ $1 == megalodon ] ; then
    mkdir -p $ssddir/megalodon

    megalodon \
        $ssddir/raw/$prefix/$prefix/fast5 \
        --overwrite \
        --guppy-server-path "/usr/bin/guppy_basecall_server" \
        --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --reference $asm \
        --outputs per_read_mods \
        --output-directory $ssddir/megalodon/${prefix}_meta \
        --write-mods-text \
        --devices "cuda:0" \
        --processes 36
fi

if [ $1 == move_meta ] ; then
    mkdir -p $datadir/megalodon/${prefix}_meta

    mv $ssddir/megalodon/${prefix}_meta/* $datadir/megalodon/${prefix}_meta/
fi


if [ $1 == megaidx ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt \
	    -o $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt.idx
fi

if [ $1 == barcode_common ] ; then
    mkdir -p $datadir/barcode/${prefix}_meta
    { time python ~/Code/yfan_meth/utils/megalodon_barcode.py \
           -m $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt \
           -i $datadir/megalodon/${prefix}_meta/per_read_modified_base_calls.txt.idx \
           -r $asm \
           -b ~/Code/yfan_nanopore/mdr/rebase/barcodes15.txt \
           -o $datadir/barcode/${prefix}_meta/${prefix}_contigs_barcodes15.txt \
           -t 36 ;} &> $datadir/barcode/${prefix}_meta/${prefix}_contigs_barcodes15_time.txt
fi
    

