#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar
ref=$datadir/reference/candida_nivariensis.fa
asm=$datadir/paperfigs/assembly_final/nivar.final.fasta

if [ $1 == final_busco ] ; then
    python ~/software/busco/scripts/run_BUSCO.py -f -i $ref -o ref -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome
    python ~/software/busco/scripts/run_BUSCO.py -f -i $asm -o asm -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome

    mkdir -p $datadir/paperfigs/busco
    mv run_asm $datadir/paperfigs/busco/
    mv run_ref $datadir/paperfigs/busco/
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


if [ $1 == transcriptome_busco ] ; then
    transcriptome=$r9dir/trinity/Trinity.fasta
    python ~/software/busco/scripts/run_BUSCO.py -f -i $transcriptome -o tran_busco -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m tran
fi

if [ $1 == augustus_busco ] ; then
    anndir=$r9dir/annotation
    
    awk '$3 == "transcript" {print $0}' $r9dir/annotation/augustus/cani_by_caal.gff > $r9dir/annotation/augustus/cani_by_caal_xscripts.gff
    bedtools getfasta -fi $anndir/nivar_fb15_bwa.fasta -bed $anndir/augustus/cani_by_caal_xscripts.gff -fo $anndir/augustus/cani_by_caal_xscripts.fasta

    python ~/software/busco/scripts/run_BUSCO.py -f -i $anndir/augustus/cani_by_caal_xscripts.fasta \
	   -o augustus_busco \
	   -l ~/software/busco/lineages/fungi_odb9 \
	   -sp candida_albicans -m tran
fi

if [ $1 == raw_busco ] ; then
    ##try with raw
    python ~/software/busco/scripts/run_BUSCO.py -f -i $r9dir/canu/nivar.contigs.fasta -o r9raw_busco -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m tran
    python ~/software/busco/scripts/run_BUSCO.py -f -i $r10dir/canu/nivar.contigs.fasta -o r10raw_busco -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m tran
fi




