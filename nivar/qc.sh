#!/bin/bash



if [ $1 == yields ] ; then
    ##get yield and trimmed yields for all data
    npr9=`bash ~/Code/utils/qc/basic_run_assess.sh $r9dir/nivar_r9.fq`
    npr9_filt=`bash ~/Code/utils/qc/basic_run_assess.sh $r9dir/nivar_r9_3k.fq`
    npr9_filt=`bash ~/Code/utils/qc/basic_run_assess.sh $r9dir/nivar_dRNA.fq`
    npr10=`bash ~/Code/utils/qc/basic_run_assess.sh $r10dir/190706_nivar_r10.fq`
    npr10_filt=`bash ~/Code/utils/qc/basic_run_assess.sh $r10dir/190706_nivar_r10_3k.fq`
    echo $npr9
    echo $npr9_filt
    echo $npr10
    echo $npr10_filt
    for i in /kyber/Data/seqlab/sp_2019/nivar_r9/data/gDNA_illumina/CANI_gDNA_R*.fastq.gz ;
    do
	yield=`zcat $i | paste - - - - | cut -f 2 | tr -d '\n' | wc -c `
	echo $i,$yield
    done
    for i in /kyber/Data/seqlab/sp_2019/nivar_r9/data/cDNA_illumina/CANI_cDNA_R*.fastq.gz ;
    do
	yield=`zcat $i | paste - - - - | cut -f 2 | tr -d '\n' | wc -c `
	echo $i,$yield
    done
    for i in /kyber/Data/seqlab/sp_2019/nivar_r9/trimmed/CANI_*_paired.fq.gz ;
    do
	yield=`zcat $i | paste - - - - | cut -f 2 | tr -d '\n' | wc -c `
	echo $i,$yield
    done
    
fi

if [ $1 == asm_stats ] ; then
    ref=$r9dir/References/candida_nivariensis.fa

    python ~/Code/utils/qc/asm_assess.py -i $ref
    python ~/Code/utils/qc/asm_assess.py -i $r9asm
    python ~/Code/utils/qc/asm_assess.py -i $r10asm
fi

if [ $1 == busco ] ; then
    r9asmdir=$r9dir/assemblies
    r10asmdir=$r10dir/assemblies

    ref=$r9dir/References/candida_nivariensis.fa

    mkdir -p $r9dir/busco
    for i in $r9asmdir/*fa* ;
    do
	prefix=`basename $i .fasta`
	python ~/software/busco/scripts/run_BUSCO.py -f -i $i -o r9_$prefix -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome
    done
    mkdir -p $r10dir/busco
    for i in $r10asmdir/*fa* ;
    do
	prefix=`basename $i .fasta`
	python ~/software/busco/scripts/run_BUSCO.py -f -i $i -o r10_$prefix -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome
    done
fi


if [ $1 == ref_busco ] ; then
    ref=$r9dir/References/candida_nivariensis.fa
    python ~/software/busco/scripts/run_BUSCO.py -f -i $ref -o ref_busco -l ~/software/busco/lineages/fungi_odb9 -sp candida_albicans -m genome
fi


if [ $1 == asmstats ] ; then
    outfile=~/Dropbox/yfan/nivar/asm/asmstats.csv
    echo asm,num,n50,longest,shortest,total > $outfile
    refstats=`python ~/Code/utils/qc/asm_assess.py -i $r9dir/References/candida_nivariensis.fa`
    echo ref,$refstats >> $outfile
    r9stats=`python ~/Code/utils/qc/asm_assess.py -i $r9dir/canu/nivar.contigs.fasta`
    echo r9,$r9stats >> $outfile
    r10stats=`python ~/Code/utils/qc/asm_assess.py -i $r10dir/canu/nivar.contigs.fasta`
    echo r10,$r10stats >> $outfile
fi


if [ $1 == mummer_corr ] ; then
    mkdir -p $r9dir/mummer_corr

    for i in $r9dir/assemblies/*fasta ;
    do
	rm ~/tmp/*
	prefix=`basename $i .fasta`
	mkdir -p $r9dir/mummer_corr/$prefix
	cp $i ~/tmp/
	cp $r9dir/canu/nivar.contigs.fasta ~/tmp/
	nucmer -p ~/tmp/$prefix ~/tmp/nivar.contigs.fasta ~/tmp/$prefix.fasta 
	mummerplot --filter --fat --png -p ~/tmp/$prefix ~/tmp/$prefix.delta
	dnadiff -p ~/tmp/$prefix ~/tmp/nivar.contigs.fasta ~/tmp/$prefix.fasta 
	cp ~/tmp/$prefix* $r9dir/mummer_corr/$prefix/
    done

    for i in $r10dir/assemblies/*fasta ;
    do
	rm ~/tmp/*
	prefix=`basename $i .fasta`
	mkdir -p $r10dir/mummer_corr/$prefix
	cp $i ~/tmp/
	cp $r10dir/canu/nivar.contigs.fasta ~/tmp/
	nucmer -p ~/tmp/$prefix ~/tmp/nivar.contigs.fasta ~/tmp/$prefix.fasta 
	mummerplot --filter --fat --png -p ~/tmp/$prefix ~/tmp/$prefix.delta
	dnadiff -p ~/tmp/$prefix ~/tmp/nivar.contigs.fasta ~/tmp/$prefix.fasta 
	cp ~/tmp/$prefix* $r10dir/mummer_corr/$prefix/
    done
fi




if [ $1 == find_enrich ] ; then
    dbxdir=~/Dropbox/yfan/nivar/error
    mkdir -p $dbxdir

    ref=$r9dir/canu/nivar.contigs.fasta
    for i in $r9dir/assemblies/*fasta ;
    do
	prefix=`basename $i .fasta`
	python2 ~/Code/utils/motif_enrich.py -s $r9dir/mummer_corr/$prefix/$prefix.snps -r $ref -m 6 -o $dbxdir/$prefix.r9.csv
    done

    
    ref=$r10dir/canu/nivar.contigs.fasta
    for i in $r10dir/assemblies/*fasta ;
    do
	prefix=`basename $i .fasta`
	python2 ~/Code/utils/motif_enrich.py -s $r10dir/mummer_corr/$prefix/$prefix.snps -r $ref -m 6 -o $dbxdir/$prefix.r10.csv
    done
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
