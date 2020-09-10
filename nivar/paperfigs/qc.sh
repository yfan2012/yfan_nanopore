#!/bin/bash

datadir=/uru/Data/Nanopore/projects/nivar

if [ $1 == final_busco ] ; then
    i=r9
    ref=$datadir/reference/candida_nivariensis.fa
    python ~/software/busco/scripts/run_BUSCO.py -f -i $datadir/pilon/${i}_pilon/nivar_${i}.pilon_bwa.6.fasta -o asm_busco -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome
    python ~/software/busco/scripts/run_BUSCO.py -f -i $ref -o ref_busco -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome

    mkdir -p $datadir/paperfigs/busco
    mv run_asm_busco $datadir/paperfigs/busco/
    mv run_ref_busco $datadir/paperfigs/busco/
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
    ref=`python ~/Code/utils/qc/asm_assess.py -i $datadir/reference/candida_nivariensis.fa`
    echo $datadir/reference/candida_nivariensis.fa,$ref
    r9=`python ~/Code/utils/qc/asm_assess.py -i $datadir/assemble/r9_assembly/nivar_r9.contigs.fasta`
    echo $datadir/assemble/r9_assembly/nivar_r9.contigs.fasta,$r9
    r10=`python ~/Code/utils/qc/asm_assess.py -i $datadir/assemble/r10_assembly/nivar_r10.contigs.fasta`
    echo $datadir/assemble/r10_assembly/nivar_r10.contigs.fasta,$r10
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

if [ $1 == extract_quals ] ; then
    python error_basequals.py -m $datadir/mummer/r9_raw_ref/nivar_r9_raw_ref.snps \
	   -b1 $datadir/align/r9/reference_r9.sorted.bam \
	   -b2 $datadir/align/r10/reference_r10.sorted.bam \
	   -o $datadir/error_quals/r9_quals.csv
    python error_basequals.py -m $datadir/mummer/r10_raw_ref/nivar_r10_raw_ref.snps \
	   -b1 $datadir/align/r10/reference_r10.sorted.bam \
	   -b2 $datadir/align/r9/reference_r9.sorted.bam \
	   -o $datadir/error_quals/r10_quals.csv
fi


