#!/bin/bash

projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/paperfigs/nohuman
prefix=200708_mdr_stool16native_nohuman


if [ $1 == assembly ] ; then
    mkdir -p $projdir/mdr/flye/$prefix

    flye \
	--nano-raw $projdir/mdr/fastqs/$prefix.fq.gz \
	-o $projdir/mdr/flye/$prefix \
	-t 36 \
	-g 100m \
	--meta
    
    mv $projdir/mdr/flye/$prefix/assembly.fasta $projdir/mdr/flye/$prefix/$prefix.assembly.fasta
fi


fq=$projdir/mdr/fastqs/$prefix.fq.gz
if [ $1 == racon_align ] ; then
    mkdir -p $projdir/mdr/racon

    minimap2 -t 36 -ax map-ont $projdir/mdr/flye/$prefix/$prefix.assembly.fasta $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $projdir/mdr/racon/$prefix.sorted.bam
    samtools index $projdir/mdr/racon/$prefix.sorted.bam
fi


if [ $1 == racon ] ; then
    mkdir -p $projdir/mdr/racon

    samtools view -@ 36 $projdir/mdr/racon/$prefix.sorted.bam > $projdir/mdr/racon/$prefix.sorted.sam
    racon -m 8 -x -6 -g -8 -w 500 -t 36 \
	  $fq \
	  $projdir/mdr/racon/$prefix.sorted.sam \
	  $projdir/mdr/flye/$prefix/$prefix.assembly.fasta > $projdir/mdr/racon/$prefix.racon.contigs.fasta
fi


if [ $1 == medaka ] ; then
    mkdir -p $projdir/mdr/medaka
    mkdir -p $projdir/mdr/medaka/$prefix

    medaka_consensus \
	-i $fq \
	-d $projdir/mdr/racon/$prefix.racon.contigs.fasta \
	-o $projdir/mdr/medaka/$prefix \
	-t 36 \
	-m r941_min_high_g360
fi

asmpolished=$projdir/mdr/medaka/$prefix/consensus.fasta
if [ $1 == align_polished ] ; then
    mkdir -p $projdir/mdr/align

    minimap2 -t 36 -ax map-ont $asmpolished $fq \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $projdir/mdr/align/$prefix.sorted.bam
    samtools index $projdir/mdr/align/$prefix.sorted.bam
fi

if [ $1 == coverage_polished ] ; then
    bedtools genomecov -d \
	     -ibam $projdir/mdr/align/$prefix.sorted.bam \
	     > $projdir/mdr/align/$prefix.sorted.cov
fi


if [ $1 == perf_aligns_only ] ; then
    ##filters for read ids with only one alignment recorded
    samtools view $projdir/mdr/align/$prefix.sorted.bam \
	| awk '{print $1}' \
	| sort \
	| uniq -u > $projdir/mdr/align/${prefix}_unique_readids.txt

    ##filters for read ids with mapq60
    samtools view $projdir/mdr/align/$prefix.sorted.bam \
	| awk '$5==60 {print $1}'\
	| sort \
	| uniq -u > $projdir/mdr/align/${prefix}_mapq60_readids.txt

    ##list of mapq60 and unqiue
    sort $projdir/mdr/align/${prefix}_unique_readids.txt $projdir/mdr/align/${prefix}_mapq60_readids.txt \
	| uniq -d >$projdir/mdr/align/${prefix}_perf_readids.txt

fi


if [ $1 == align_perf ] ; then
    seqtk subseq $fq $projdir/mdr/align/${prefix}_perf_readids.txt > $projdir/mdr/fastqs/$prefix.perf.fq
    pigz -p 12 $projdir/mdr/fastqs/$prefix.perf.fq

    minimap2 -t 36 -ax map-ont $asmpolished $projdir/mdr/fastqs/$prefix.perf.fq.gz \
	| samtools view -@ 36 -b \
	| samtools sort -@ 36 -o $projdir/mdr/align/$prefix.perf.sorted.bam
    samtools index $projdir/mdr/align/$prefix.perf.sorted.bam


    minimap2 -t 36 -ax map-ont $asmpolished $projdir/mdr/fastqs/$prefix.perf.fq.gz \
	     > $projdir/mdr/align/$prefix.perf.paf

    bedtools genomecov -d \
	     -ibam $projdir/mdr/align/$prefix.perf.sorted.bam \
	     > $projdir/mdr/align/$prefix.perf.sorted.cov
fi


ssddir=~/data/mdr/mdr
if [ $1 == find_perf_fast5 ] ; then
    fast5_subset \
	-i $ssddir/raw/$prefix/fast5 \
	-s $ssddir/raw/${prefix}_perf \
	-l $projdir/mdr/align/${prefix}_perf_readids.txt \
	--recursive
fi


if [ $1 == megalodon_perf ] ; then
    ##run megalodon with the asm and reads that only align uniquely with perf mapq
    megalodon \
	$ssddir/raw/${prefix}_perf \
	--overwrite \
	--suppress-progress-bars \
	--verbose-read-progress 0 \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	--reference $asmpolished \
	--outputs per_read_mods mod_basecalls \
	--output-directory $ssddir/megalodon/${prefix}_perf \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36
    rm $ssddir/megalodon/${prefix}_perf/per_read_modified_base_calls.db
    
    mkdir -p $projdir/mdr/megalodon/${prefix}_perf
    mv $ssddir/megalodon/${prefix}_perf/* $projdir/mdr/megalodon/${prefix}_perf/
fi


if [ $1 == megaidx ] ; then
    python3 ~/Code/yfan_meth/utils/megalodon_mod_basecalls_idx.py \
	    -i $projdir/mdr/megalodon/${prefix}_perf/per_read_modified_base_calls.txt \
	    -o $projdir/mdr/megalodon/${prefix}_perf/per_read_modified_base_calls.txt.idx
fi


barcodes3=~/Code/yfan_nanopore/mdr/paperfigs/contig_level/clin_barcodes3.txt
if [ $1 == methprobs ] ; then
    mkdir -p $projdir/paperfigs/nohuman
    
    ##filter out barcode motif related positions in per read meth calls
    python3 ~/Code/yfan_meth/utils/megalodon_extract_barcode_methprobs.py \
	   -r $asmpolished \
	   -b $barcodes3 \
	   -m $projdir/mdr/megalodon/${prefix}_perf/per_read_modified_base_calls.txt \
	   -i $projdir/mdr/megalodon/${prefix}_perf/per_read_modified_base_calls.txt.idx \
	   -o $datadir/methprobs_nohuman_perf.csv \
	   -t 12

fi


if [ $1 == methcalls ] ; then
    ##assign meth/unmeth based on given thresholds
    python3 ~/Code/yfan_nanopore/mdr/zymo/contig_agg/filter_motif_calls.py \
	   -i $datadir/methprobs_nohuman_perf.csv \
	   -o $datadir/methcalls_nohuman_perf.csv \
	   -m .8 \
	   -u .8
fi


if [ $1 == mummer_hiC ] ; then
    ##mummer clin sample contigs with bins from phase folks
    mkdir -p $datadir/clin_mummer
    
    nucmer -p $datadir/clin_mummer/asm_hiC $asmpolished $projdir/mdr/hiC/200708_mdr_stool16native.hiC.fasta
    mummerplot --filter --fat --postscript -p $datadir/clin_mummer/asm_hiC $datadir/clin_mummer/asm_hiC.delta
    mummerplot --filter --fat --png -p $datadir/clin_mummer/asm_hiC $datadir/clin_mummer/asm_hiC.delta
    dnadiff -p $datadir/clin_mummer/asm_hiC $asmpolished $projdir/mdr/hiC/200708_mdr_stool16native.hiC.fasta
fi


if [ $1 == amr ] ; then
    mkdir -p $projdir/mdr/amr

    for i in ecoh card ncbi resfinder plasmidfinder vfdb ecoli_vf megares argannot ;
    do
	abricate \
	    --threads 36 \
	    --db $i \
	    $asmpolished > $projdir/mdr/amr/$prefix.$i.tsv
    done
fi


if [ $1 == asm_assess ] ; then
    python3 ~/Code/utils/qc/asm_assess.py \
	    -i $asmpolished \
	    -p mdr > ~/gdrive/mdr/paperfigs/figs/asmstats.csv
fi


catdir=/atium/Data/ref/CATdb
if [ $1 == cattigs ] ; then
    mkdir -p $datadir/contig_id/CAT

    CAT contigs \
	-c $asmpolished \
	-d $catdir/db \
	-t $catdir/tax \
	-o $datadir/contig_id/CAT/$prefix.CAT \
	--force \
	--sensitive

    CAT add_names \
	-i $datadir/contig_id/CAT/$prefix.CAT.contig2classification.txt \
	-t $catdir/tax \
	-o $datadir/contig_id/CAT/$prefix.CAT.names.txt

    CAT add_names \
	-i $datadir/contig_id/CAT/$prefix.CAT.contig2classification.txt \
	-t $catdir/tax \
	-o $datadir/contig_id/CAT/$prefix.CAT.names_official.txt \
	--only_official
    
    CAT summarise \
	-c $asmpolished \
	-i $datadir/contig_id/CAT/$prefix.CAT.names_official.txt \
	-o $datadir/contig_id/CAT/$prefix.CAT.summary.txt
fi


if [ $1 == kraken ] ; then
    dbdir=/mithril/Data/Nanopore/ref/kraken2
    
    krakendir=$datadir/contig_id/kraken_standard
    mkdir -p $krakendir
    kraken2 \
	--db $dbdir/standard \
	--threads 36 \
	--classified-out $krakendir/$prefix.class.txt \
	--unclassified-out $krakendir/$prefix.unclass.txt \
	--output $krakendir/$prefix.out.txt \
	--report $krakendir/$prefix.report.txt \
	--use-names \
	$asmpolished

    krakenont=$datadir/contig_id/kraken_standard_ont
    mkdir -p $krakenont
    kraken2 \
	--db $dbdir/standard_ont \
	--threads 36 \
	--classified-out $krakenont/$prefix.class.txt \
	--unclassified-out $krakenont/$prefix.unclass.txt \
	--output $krakenont/$prefix.out.txt \
	--report $krakenont/$prefix.report.txt \
	--use-names \
	$asmpolished

fi
