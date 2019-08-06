#!/bin/bash

datadir=~/Dropbox/yfan/nina_fungus/assemblies
tmp=~/tmp/assemblies

if [ $1 == mum_assemblers ] ; then
    ##mummer assemblies against each other. Try to ditch wtdbg2 junk contigs in some way
    mkdir -p ~/mummer_assemblers

    ##clean up contig names in fastas and copy to a place where mummer can find it
    sed -i -e 's/_pilon//g' $datadir/*pilon/*fasta
    rm -r ~/tmp
    mkdir ~/tmp
    cp -r $datadir ~/tmp/
    
    for i in st31 st90853 ;
    do
	##mummer is being stupid about file paths. so write results to home first, then move the folder when everything is done
	nucmer -p ~/mummer_assemblers/$i $tmp/wtdbg2_pilon/$i*.fasta $tmp/canu_pilon/$i*20.fasta

	mummerplot --filter --fat --png -p ~/mummer_assemblers/$i.layout ~/mummer_assemblers/$i.delta -R $tmp/wtdbg2_pilon/$i*.fasta -Q $tmp/canu_pilon/$i*20.fasta
	mummerplot --filter --fat --png -p ~/mummer_assemblers/$i ~/mummer_assemblers/$i.delta

	dnadiff -p ~/mummer_assemblers/$i $tmp/wtdbg2_pilon/$i*.fasta $tmp/canu_pilon/$i*20.fasta
    done
    rm -rf $datadir/mummer_assemblers
    mv ~/mummer_assemblers $datadir/
fi

    
if [ $1 == align_np ] ; then
    ##align nanopore reads to check misassembly
    datadir=/scratch/groups/mschatz1/cpowgs/fungus/181108_nina_v2

    mkdir -p $datadir/align

    for i in st31 st90853 ;
    do
	minimap2 -a -x map-ont -t 36 $datadir/$i/canu_pilon_trimmed_$i/${i}_canu.pilon.19.fasta $datadir/fastqs/${i}_bothruns_over3kb.fastq | samtools view -b | samtools sort -o $datadir/align/${i}_canu.pilon.19.sorted.bam -T $datadir/align/reads.tmp -
	samtools index $datadir/align/${i}_canu.pilon.19.sorted.bam
    done
    
    minimap2 -a -x map-ont -t 36 $datadir/st31/pilon_trimmed_st31/st31_wtdbg2.pilon.20.fasta $datadir/fastqs/st31_bothruns_over3kb.fastq | samtools view -b | samtools sort -o $datadir/align/st31_wtdbg2.pilon.20.sorted.bam -T $datadir/align/reads.tmp -
    samtools index $datadir/align/st31_wtdbg2.pilon.20.sorted.bam
    
    minimap2 -a -x map-ont -t 36 $datadir/st90853/pilon_trimmed_st90853/st90853_wtdbg2.pilon.25.fasta $datadir/fastqs/st90853_bothruns_over3kb.fastq | samtools view -b | samtools sort -o $datadir/align/st90853_wtdbg2.pilon.25.sorted.bam -T $datadir/align/reads.tmp -
    samtools index $datadir/align/st90853_wtdbg2.pilon.25.sorted.bam
fi
