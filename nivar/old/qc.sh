#!/bin/bash

r9dir=/kyber/Data/seqlab/sp_2019/nivar_r9
r10dir=/kyber/Data/Nanopore/Analysis/190706_nivar_r10

r9asm=$r9dir/canu/nivar.contigs.fasta
r10asm=$r10dir/canu/nivar.contigs.fasta

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

if [ $1 == asms ] ; then
    ref=$r9dir/References/candida_nivariensis.fa

    python ~/Code/utils/qc/asm_assess.py -i $ref
    python ~/Code/utils/qc/asm_assess.py -i $r9asm
    python ~/Code/utils/qc/asm_assess.py -i $r10asm
fi


if [ $1 == ref_mum ] ; then
    mkdir -p $r9dir/mummer
    mkdir -p $r10dir/mummer

    cp $r9asm ~/tmp/
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r9 ~/tmp/nivar.contigs.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r9 ~/tmp/nivar_r9.delta
    dnadiff -p ~/tmp/nivar_r9 ~/tmp/nivar.contigs.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r9* $r9dir/mummer/

    
    cp $r10asm ~/tmp/
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r10 ~/tmp/nivar.contigs.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r10 ~/tmp/nivar_r10.delta
    dnadiff -p ~/tmp/nivar_r10 ~/tmp/nivar.contigs.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r10* $r10dir/mummer/
fi

if [ $1 == clean_pilon ] ; then
    ##clean up pilon outputs
    sed -i -e 's/_pilon//g' $r9dir/pilon_bwa/nivar*15*fasta
    sed -i -e 's/>/>pilon_10_/g' $r9dir/pilon_bwa/nivar*15*fasta
    sed -i -e 's/_pilon//g' $r10dir/pilon_bwa/nivar*15*fasta
    sed -i -e 's/>/>pilon_10_/g' $r10dir/pilon_bwa/nivar*15*fasta
fi



if [ $1 == error_mum_pilon ] ; then
    rm -r ~/tmp/*
    mkdir -p $r9dir/mummer/pilon
    cp $r9dir/pilon_bwa/nivar.pilon_bwa.15.fasta ~/tmp
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r9_pilon_bwa ~/tmp/nivar.pilon_bwa.15.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r9_pilon_bwa ~/tmp/nivar_r9_pilon_bwa.delta
    dnadiff -p ~/tmp/nivar_r9_pilon_bwa ~/tmp/nivar.pilon_bwa.15.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r9_pilon* $r9dir/mummer/pilon/
    
    rm -r ~/tmp/*
    mkdir -p $r10dir/mummer/pilon
    cp $r10dir/pilon_bwa/nivar.pilon_bwa.15.fasta ~/tmp
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r10_pilon_bwa ~/tmp/nivar.pilon_bwa.15.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r10_pilon_bwa ~/tmp/nivar_r10_pilon_bwa.delta
    dnadiff -p ~/tmp/nivar_r10_pilon_bwa ~/tmp/nivar.pilon_bwa.15.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r10_pilon* $r10dir/mummer/pilon/
fi

if [ $1 == error_mum_racon ] ; then
    rm -r ~/tmp/*
    mkdir -p $r9dir/mummer/racon
    cp $r9dir/racon/nivar_bwa.racon.15.fasta ~/tmp/
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r9_racon_bwa ~/tmp/nivar_bwa.racon.15.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r9_racon_bwa ~/tmp/nivar_r9_racon_bwa.delta
    dnadiff -p ~/tmp/nivar_r9_racon_bwa ~/tmp/nivar_bwa.racon.15.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r9_racon_bwa* $r9dir/mummer/racon/


    rm -r ~/tmp/*
    mkdir -p $r10dir/mummer/racon
    cp $r10dir/racon/nivar_bwa.racon.15.fasta ~/tmp/
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r10_racon_bwa ~/tmp/nivar_bwa.racon.15.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r10_racon_bwa ~/tmp/nivar_r10_racon_bwa.delta
    dnadiff -p ~/tmp/nivar_r10_racon_bwa ~/tmp/nivar_bwa.racon.15.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r10_racon_bwa* $r10dir/mummer/racon/
fi


if [ $1 == error_mum_freebayes ] ; then
    rm -r ~/tmp/*
    mkdir -p $r9dir/mummer/freebayes
    cp $r9dir/freebayes_bwa/nivar_fb15_bwa.fasta ~/tmp/
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r9_fb_bwa ~/tmp/nivar_fb15_bwa.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r9_fb_bwa ~/tmp/nivar_r9_fb_bwa.delta
    dnadiff -p ~/tmp/nivar_r9_fb_bwa ~/tmp/nivar_fb15_bwa.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r9_fb_bwa* $r9dir/mummer/freebayes/


    rm -r ~/tmp/*
    mkdir -p $r10dir/mummer/freebayes
    cp $r10dir/freebayes_bwa/nivar_fb15_bwa.fasta ~/tmp/
    cp $r9dir/References/candida_nivariensis.fa ~/tmp/
    nucmer -p ~/tmp/nivar_r10_fb_bwa ~/tmp/nivar_fb15_bwa.fasta ~/tmp/candida_nivariensis.fa
    mummerplot --filter --fat --png -p ~/tmp/nivar_r10_fb_bwa ~/tmp/nivar_r10_fb_bwa.delta
    dnadiff -p ~/tmp/nivar_r10_fb_bwa ~/tmp/nivar_fb15_bwa.fasta ~/tmp/candida_nivariensis.fa
    cp ~/tmp/nivar_r10_fb_bwa* $r10dir/mummer/freebayes/
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
