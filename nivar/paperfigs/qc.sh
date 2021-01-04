#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar

ref=$datadir/reference/candida_nivariensis.fa
gla=$datadir/reference/candida_glabrata.fa
glagff=$datadir/reference/candida_glabrata.gff
glatransfa=$datadir/reference/candida_glabrata.transcriptome.fa
cer=$datadir/reference/saccharomyces_cerevisiae.fa
cergff=$datadir/reference/saccharomyces_cerevisiae.gff
certransfa=$datadir/reference/saccharomyces_cerevisiae.transcriptome.fa
alb=$datadir/reference/candida_albicans.fa
albgff=$datadir/reference/candida_albicans.gff
albtransfa=$datadir/reference/candida_albicans.transcriptome.fa

asm=$datadir/paperfigs/assembly_final/nivar.final.fasta
asmgff=$datadir/paperfigs/annotation_final/nivar.final.gff
asmtransfa=$datadir/paperfigs/annotation_final/nivar.final.transcriptome.fasta

bragff=$datadir/paperfigs/annotation/braker/braker.gff3
bratransfa=$datadir/paperfigs/annotation/braker/braker.fa

fqraw=$datadir/r9/r9.fq
fq=$datadir/r9/r9_3kb.fq
illfwd=$datadir/illumina/gDNA_trimmed/nivar_gDNA_fwd_paired.fq.gz
illrev=$datadir/illumina/gDNA_trimmed/nivar_gDNA_rev_paired.fq.gz

if [ $1 == final_busco ] ; then
    busco \
	-m genome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $asm \
	-o asm \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f

    busco \
	-m genome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $ref \
	-o ref \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f 
fi

if [ $1 == ref_busco ] ; then
    busco \
	-m genome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $gla \
	-o gla \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f

    busco \
	-m genome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $cer \
	-o cer \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f
    
    busco \
	-m genome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $alb \
	-o alb \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f 
	
fi

if [ $1 == ref_busco_transcriptome ] ; then
    #combined
    bedtools getfasta \
	     -fi $gla \
	     -bed $glagff \
	     -fo $glatransfa

    busco \
	-m transcriptome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $glatransfa \
	-o transgla \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f
    
    bedtools getfasta \
	     -fi $cer \
	     -bed $cergff \
	     -fo $certransfa

    busco \
	-m transcriptome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $certransfa \
	-o transcer \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f
    
    bedtools getfasta \
	     -fi $alb \
	     -bed $albgff \
	     -fo $albtransfa

    busco \
	-m transcriptome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $albtransfa \
	-o transalb \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f
fi

if [ $1 == busco_transcriptome ] ; then
    #braker
    bedtools getfasta \
	     -fi $asm \
	     -bed $bragff \
	     -fo $bratransfa
    
    busco \
	-m transcriptome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $bratransfa \
	-o transasm_brakeronly \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f
fi

if [ $1 == busco_transcriptome_final ] ; then
    bedtools getfasta \
	     -fi $asm \
	     -bed $asmgff \
	     -fo $asmtransfa
    busco \
	-m transcriptome \
	-l saccharomycetes_odb10 \
	--augustus_species saccharomyces_cerevisiae_S288C \
	-i $asmtransfa \
	-o transasm \
	--out_path $datadir/paperfigs/busco \
	-c 36 \
	-f
fi



if [ $1 == yields ] ; then
    
    for i in $datadir/illumina/gDNA/nivar_gDNA_R*.fastq.gz ;
    do
	newfile=`basename $i .gz`
	zcat $i > $datadir/illumina/gDNA/$newfile
	bash ~/Code/utils/qc/basic_run_assess.sh $datadir/illumina/gDNA/$newfile
    done
    
    for i in $datadir/illumina/gDNA_trimmed/nivar_gDNA_*_paired.fq.gz ;
    do
	newfile=`basename $i .gz`
	zcat $i > $datadir/illumina/gDNA_trimmed/$newfile
	bash ~/Code/utils/qc/basic_run_assess.sh $datadir/illumina/gDNA_trimmed/$newfile
    done

    for i in $datadir/illumina/cDNA/nivar_cDNA_R*.fastq.gz ;
    do
	newfile=`basename $i .gz`
	zcat $i > $datadir/illumina/cDNA/$newfile
	bash ~/Code/utils/qc/basic_run_assess.sh $datadir/illumina/cDNA/$newfile
    done
    
    for i in $datadir/illumina/cDNA_trimmed/nivar_cDNA_*_paired.fq.gz ;
    do
	newfile=`basename $i .gz`
	zcat $i > $datadir/illumina/cDNA_trimmed/$newfile
	bash ~/Code/utils/qc/basic_run_assess.sh $datadir/illumina/cDNA_trimmed/$newfile
    done

    ##get yield and trimmed yields for all data
    npr9=`bash ~/Code/utils/qc/basic_run_assess.sh $datadir/r9/r9.fq`
    npr9_filt=`bash ~/Code/utils/qc/basic_run_assess.sh $datadir/r9/r9_3kb.fq`
    drna=`bash ~/Code/utils/qc/basic_run_assess.sh $datadir/dRNA/dRNA.fq`
    npr10=`bash ~/Code/utils/qc/basic_run_assess.sh $datadir/r10/r10.fq`
    npr10_filt=`bash ~/Code/utils/qc/basic_run_assess.sh $datadir/r10/r10_3kb.fq`

    echo $npr9
    echo $npr9_filt
    echo $npr10
    echo $npr10_filt
    echo $drna

fi

if [ $1 == asm_stats ] ; then
    ##prints asm stats for reference and final
    refstats=`python ~/Code/utils/qc/asm_assess.py -i $ref`
    echo $ref,$refstats
    asmstats=`python ~/Code/utils/qc/asm_assess.py -i $asm`
    echo $asm,$asmstats
fi


if [ $1 == mummer_final ] ; then
    mkdir -p ~/tmp/mummer
    mkdir -p ~/tmp/mummer/ref
    mkdir -p ~/tmp/mummer/cer
    mkdir -p ~/tmp/mummer/gla
    
    cp $asm ~/tmp/mummer/
	
    cp $ref ~/tmp/mummer/ref
    nucmer -p ~/tmp/mummer/ref/ref_nivar.final ~/tmp/mummer/nivar.final.fasta ~/tmp/mummer/ref/candida_nivariensis.fa 
    mummerplot --filter --fat --postscript -p ~/tmp/mummer/ref/ref_nivar.final ~/tmp/mummer/ref/ref_nivar.final.delta
    mummerplot --filter --fat --png -p ~/tmp/mummer/ref/ref_nivar.final ~/tmp/mummer/ref/ref_nivar.final.delta
    dnadiff -p ~/tmp/mummer/ref/ref_nivar.final ~/tmp/mummer/nivar.final.fasta ~/tmp/mummer/ref/candida_nivariensis.fa 

    cp $gla ~/tmp/mummer/gla
    nucmer -p ~/tmp/mummer/gla/gla_nivar.final ~/tmp/mummer/nivar.final.fasta ~/tmp/mummer/gla/candida_glabrata.fa 
    mummerplot --filter --fat --postscript -p ~/tmp/mummer/gla/gla_nivar.final ~/tmp/mummer/gla/gla_nivar.final.delta
    mummerplot --filter --fat --png -p ~/tmp/mummer/gla/gla_nivar.final ~/tmp/mummer/gla/gla_nivar.final.delta
    dnadiff -p ~/tmp/mummer/gla/gla_nivar.final ~/tmp/mummer/nivar.final.fasta ~/tmp/mummer/gla/candida_glabrata.fa 
    
    cp $cer ~/tmp/mummer/cer
    nucmer -p ~/tmp/mummer/cer/cer_nivar.final ~/tmp/mummer/nivar.final.fasta ~/tmp/mummer/cer/saccharomyces_cerevisiae.fa 
    mummerplot --filter --fat --postscript -p ~/tmp/mummer/cer/cer_nivar.final ~/tmp/mummer/cer/cer_nivar.final.delta
    mummerplot --filter --fat --png -p ~/tmp/mummer/cer/cer_nivar.final ~/tmp/mummer/cer/cer_nivar.final.delta
    dnadiff -p ~/tmp/mummer/cer/cer_nivar.final ~/tmp/mummer/nivar.final.fasta ~/tmp/mummer/cer/saccharomyces_cerevisiae.fa

fi
    
if [ $1 == coverage ] ; then
    mkdir -p $datadir/paperfigs/cov
    
    minimap2 -ax map-ont -t 36 $asm $fq |
	samtools view -@ 36 -bS |
	samtools sort -@ 36 -o $datadir/paperfigs/align/nivar.final_nanopore.sorted.bam
    samtools index $datadir/paperfigs/align/nivar.final_nanopore.sorted.bam
    bedtools genomecov -d -ibam $datadir/paperfigs/align/nivar.final_nanopore.sorted.bam > $datadir/paperfigs/cov/nivar.final_nanopore.cov
    
    bwa index $asm
    bwa mem -t 36 $asm $illfwd $illrev |
	samtools view -@ 36 -bS - |
	samtools sort -@ 36 -o $datadir/paperfigs/align/nivar.final_illumina.sorted.bam
    samtools index $datadir/paperfigs/align/nivar.final_illumina.sorted.bam
    bedtools genomecov -d -ibam $datadir/paperfigs/align/nivar.final_illumina.sorted.bam > $datadir/paperfigs/cov/nivar.final_illumina.cov
    
fi

if [ $1 == coverage_raw ] ; then
    mkdir -p $datadir/paperfigs/cov
    
    minimap2 -ax map-ont -t 36 $datadir/paperfigs/assembly/nivar.contigs.fasta $fq |
	samtools view -@ 36 -bS |
	samtools sort -@ 36 -o $datadir/paperfigs/align/nivar.contigs_nanopore.sorted.bam
    samtools index $datadir/paperfigs/align/nivar.contigs_nanopore.sorted.bam
    bedtools genomecov -d -ibam $datadir/paperfigs/align/nivar.contigs_nanopore.sorted.bam > $datadir/paperfigs/cov/nivar.contigs_nanopore.cov
    
fi

