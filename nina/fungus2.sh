#!/bin/bash

srcdir=~/Code/utils/marcc
datadir=/scratch/groups/mschatz1/cpowgs/fungus



if [ $1 == untar_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2
    mkdir -p $datadir/181108_nina_v2/raw
    mkdir -p $datadir/181108_nina_v2/batch_logs
    ##sbatch --output=$datadir/181108_nina_v2/batch_logs/untar.out --job-name=ut_nina $srcdir/untar.scr $datadir/181108_nina_v2.tar.gz $datadir/181108_nina_v2
    ##shitty lustre file system caused untar to time out
    bash $srcdir/untar.scr $datadir/181108_nina_v2.tar.gz $datadir/181108_nina_v2
fi
    


if [ $1 == call_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/called
    mkdir -p $datadir/181108_nina_v2/call_logs
    mkdir -p $datadir/181108_nina_v2/call_done
    sbatch --array=0-1901 --job-name=nina_call --output=$datadir/181108_nina_v2/call_logs/nina_call.%A_%a.out $srcdir/bc_call_LSK109.scr $datadir/181108_nina_v2
fi



if [ $1 == fastq_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/fastqs/
    cat $datadir/181108_nina_v2/called/*/workspace/pass/barcode12/*.fastq > $datadir/181108_nina_v2/fastqs/st90853.fastq
    cat $datadir/181108_nina_v2/called/*/workspace/pass/barcode11/*.fastq > $datadir/181108_nina_v2/fastqs/st31.fastq
fi



if [ $1 == assemble_90853_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/st90853_assembly
    canu \
	-p st90853 -d $datadir/181108_nina_v2/st90853_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=39m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/181108_nina_v2/fastqs/st90853_bothruns_over3kb.fastq
fi



if [ $1 == assemble_31_v2 ] ; then
    mkdir -p $datadir/181108_nina_v2/st31_assembly
    canu \
	-p st31 -d $datadir/181108_nina_v2/st31_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=39m \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/181108_nina_v2/fastqs/st31_bothruns_over3kb.fastq
fi


if [ $1 == merge_fq ] ; then
    cat $datadir/181108_nina_v2/fastqs/st90853.fastq $datadir/180827_nina_fungus2/fastqs/st90853.fastq > $datadir/181108_nina_v2/fastqs/st90853_bothruns.fastq
    cat $datadir/181108_nina_v2/fastqs/st31.fastq $datadir/180827_nina_fungus2/fastqs/st31.fastq > $datadir/181108_nina_v2/fastqs/st31_bothruns.fastq
fi
    
if [ $1 == long_fq ] ; then
    python ~/Code/utils/fastq_long.py -i $datadir/181108_nina_v2/fastqs/st31_bothruns.fastq -o $datadir/181108_nina_v2/fastqs/st31_bothruns_over3kb.fastq -l 3000
    python ~/Code/utils/fastq_long.py -i $datadir/181108_nina_v2/fastqs/st90853_bothruns.fastq -o $datadir/181108_nina_v2/fastqs/st90853_bothruns_over3kb.fastq -l 3000
fi

if [ $1 == assemble_90853_wtdbg2 ] ; then
    mkdir -p $datadir/181108_nina_v2/st90853_wtdbg2
    wtdbg2 -t 32 -i $datadir/181108_nina_v2/fastqs/st90853_bothruns_over3kb.fastq -fo $datadir/181108_nina_v2/st90853_wtdbg2/st90853_wtdbg2
    wtpoa-cns -t 32 -i $datadir/181108_nina_v2/wtdbg2/st90853_wtdbg2.ctg.lay.gz -fo $datadir/181108_nina_v2/st90853_wtdbg2/st90853.wtdbg2.contigs.fasta
fi

if [ $1 == assemble_31_wtdbg2 ] ; then
    mkdir -p $datadir/181108_nina_v2/st31_wtdbg2
    wtdbg2 -t 32 -i $datadir/181108_nina_v2/fastqs/st31_bothruns_over3kb.fastq -fo $datadir/181108_nina_v2/st31_wtdbg2/st31_wtdbg2
    wtpoa-cns -t 32 -i $datadir/181108_nina_v2/st31_wtdbg2/st31_wtdbg2.ctg.lay.gz -fo $datadir/181108_nina_v2/st31_wtdbg2/st31.wtdbg2.contigs.fasta
fi


if [ $1 == pilon ] ; then
    ##sed -i -e 's/ /_/g' $datadir/181108_nina_v2/st31_wtdbg2/st31.wtdbg2.contigs.fasta
    ##sed -i -e 's/ /_/g' $datadir/181108_nina_v2/st90853_wtdbg2/st90853.wtdbg2.contigs.fasta
    
    mkdir -p $datadir/181108_nina_v2/pilon_st90853
    mkdir -p $datadir/181108_nina_v2/pilon_st31

    ##cp /work-zfs/mschatz1/cpowgs/fungus/illumina/st31* $datadir/181108_nina_v2/pilon_st31/
    ##cp /work-zfs/mschatz1/cpowgs/fungus/illumina/st90853* $datadir/181108_nina_v2/pilon_st90853/
    
    sbatch --output=$datadir/181108_nina_v2/batch_logs/st31_wtdbg2.out --job-name=st31_wtdbg2 ./pilon.scr $datadir/181108_nina_v2/pilon_st31 $datadir/181108_nina_v2/st31_wtdbg2/st31.wtdbg2.contigs.fasta st31 wtdbg2
    sbatch --output=$datadir/181108_nina_v2/batch_logs/st90853_wtdbg2.out --job-name=st90853_wtdbg2 ./pilon.scr $datadir/181108_nina_v2/pilon_st90853 $datadir/181108_nina_v2/st90853_wtdbg2/st90853.wtdbg2.contigs.fasta st90853 wtdbg2
fi

if [ $1 == pilon2 ] ; then
    ##sed -i -e 's/ /_/g' $datadir/181108_nina_v2/st31_wtdbg2/st31.wtdbg2.contigs.fasta
    ##sed -i -e 's/ /_/g' $datadir/181108_nina_v2/st90853_wtdbg2/st90853.wtdbg2.contigs.fasta
    
    mkdir -p $datadir/181108_nina_v2/pilon_st90853
    mkdir -p $datadir/181108_nina_v2/pilon_st31

    ##cp /work-zfs/mschatz1/cpowgs/fungus/illumina/st31* $datadir/181108_nina_v2/pilon_st31/
    ##cp /work-zfs/mschatz1/cpowgs/fungus/illumina/st90853* $datadir/181108_nina_v2/pilon_st90853/
    
    sbatch --output=$datadir/181108_nina_v2/batch_logs/st31_wtdbg2.out --job-name=st31_wtdbg2 ./pilon2.scr $datadir/181108_nina_v2/pilon_st31 $datadir/181108_nina_v2/st31_wtdbg2/st31.wtdbg2.contigs.fasta st31 wtdbg2
    sbatch --output=$datadir/181108_nina_v2/batch_logs/st90853_wtdbg2.out --job-name=st90853_wtdbg2 ./pilon2.scr $datadir/181108_nina_v2/pilon_st90853 $datadir/181108_nina_v2/st90853_wtdbg2/st90853.wtdbg2.contigs.fasta st90853 wtdbg2
fi
