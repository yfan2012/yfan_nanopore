#!/bin/bash

datadir=/kyber/Data/seqlab/sp_2019/fungus_asm

if [ $1 == augustus ] ; then
    anndir=$datadir/annotation
    mkdir -p $anndir/augustus
    
    augustus --species=candida_albicans $anndir/candida_nivariensis_canu.pilon.6.fasta > $anndir/augustus/cani_by_caal.gff
fi

if [ $1 == medusa ] ; then
    scadir=$datadir/medusa
    mkdir -p $scadir

    cp -r ~/software/medusa/medusa_scripts ./
    
    java -jar ~/software/medusa/medusa.jar -f $datadir/References -i $datadir/annotation/candida_nivariensis_canu.pilon.6.fasta -v -o $scadir/cani_scaffold.fasta
fi

if [ $1 == mito_blast ] ; then
    ##downloaded http://mitofun.biol.uoa.gr/fasta/sCDS.fasta.zip to find mito seq in scaffolds
    makeblastdb -in $datadir/medusa/cani_scaffold.fasta -out $datadir/blast/cani_scaffold_db -dbtype nucl
    blastn -query $datadir/References/mito/all_genes.fa -db $datadir/blast/cani_scaffold_db -outfmt 7 -out $datadir/blast/mitohits.tsv
fi
    
if [ $1 == gather_mito ] ; then
    grep Scaffold $datadir/blast/mitohits.tsv > $datadir/blast/mitohits_only.tsv
fi


if [ $1 == busco ] ; then
    mkdir -p $datadir/busco
    mkdir -p $datadir/busco_ref
    cd $datadir/busco
    run_BUSCO.py -c 12 -f -i $datadir/medusa/cani_scaffold.fasta -o scaffold_fundb -l ~/software/busco/lineages/fungi_odb9 -m geno
    cd $datadir/busco_ref
    run_BUSCO.py -c 12 -f -i $datadir/References/candida_nivariensis.fa -o scaffold_fundb -l ~/software/busco/lineages/fungi_odb9 -m geno
fi

if [ $1 == trinity ] ; then
    mkdir -p $datadir/trinity
    Trinity --seqType fq \
	    --left $datadir/data/trimmed/CANI_cDNA_forward_paired.fq.gz \
	    --right $datadir/data/trimmed/CANI_cDNA_reverse_paired.fq.gz \
	    --CPU 32 --max_memory 100G --output $datadir/trinity
fi

if [ $1 == maker_prep ] ; then
    mkdir -p $datadir/maker
    cd $datadir/maker
    ##make control files for maker
    maker -CTL
fi


if [ $1 == maker_contigs ] ; then
    cd $datadir/maker
    maker -c 36 -f -g $datadir/canu_pilon/candida_nivariensis_canu.pilon.6.fasta \
	  -base contigs \
	  $datadir/maker/maker_opts.ctl $datadir/maker/maker_bopts.ctl $datadir/maker/maker_exe.ctl
fi

##run with this example https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
if [ $1 == repeats ] ; then
    mkdir -p $datadir/repeats
    cd $datadir/repeats
    BuildDatabase -name cani -engine ncbi $datadir/medusa/cani_scaffold.fasta
    RepeatModeler -pa 36 -engine ncbi -database cani 2>&1 | tee repeatmodeler.log
fi   
   

if [ $1 == maker_scaffold ] ; then
    cd $datadir/maker
    maker -c 36 -f -g $datadir/medusa/cani_scaffold.fasta \
	  -base scaffold \
	  $datadir/maker/maker_opts.ctl $datadir/maker/maker_bopts.ctl $datadir/maker/maker_exe.ctl
fi


if [ $1 == trinity_dRNA ] ; then
    datadir=/kyber/Data/seqlab/sp_2019/nivar_r9
    mkdir -p $datadir/trinity_dRNA
    Trinity --seqType fq \
	    --single $datadir/nivar_dRNA.fq \
	    --CPU 32 --max_memory 100G --output $datadir/trinity_dRNA
fi


if [ $1 == map_drna ] ; then
    datadir=/kyber/Data/seqlab/sp_2019/nivar_r9
    mkdir -p $datadir/rna_align
    minimap2 -a -x splice -uf -k14 -t 36 $datadir/freebayes_bwa/nivar_fb15_bwa.fasta $datadir/nivar_dRNA.fq | \
	samtools view -@ 36 -b | \
	samtools sort -@ 36 -o $datadir/rna_align/nivar_fb15_dRNA.sorted.bam -T $datadir/rna_align/reads.tmp -
    samtools index $datadir/rna_align/nivar_fb15_dRNA.sorted.bam
fi
    
if [ $1 == augustus_fa ] ; then
    datadir=/kyber/Data/seqlab/sp_2019/nivar_r9
    anndir=$datadir/annotation
    cp $datadir/freebayes_bwa/nivar_fb15_bwa.fasta $anndir/
    augustus --species=candida_albicans $anndir/nivar_fb15_bwa.fasta > $anndir/augustus/cani_by_caal.gff
    bedtools getfasta -fi $anndir/nivar_fb15_bwa.fasta -bed $anndir/augustus/cani_by_caal.gff -fo $anndir/augustus/cani_by_caal.fasta
fi


if [ $1 == map_cdna ] ; then
    datadir=/kyber/Data/seqlab/sp_2019/nivar_r9

    hisat2-build $datadir/freebayes_bwa/nivar_fb15_bwa.fasta $datadir/freebayes_bwa/nivar_fb15_bwa
    hisat2 -p 36 -q -x $datadir/freebayes_bwa/nivar_fb15_bwa \
	   -1 $datadir/trimmed/CANI_cDNA_forward_paired.fq.gz \
	   -2 $datadir/trimmed/CANI_cDNA_reverse_paired.fq.gz | \
	samtools view -@ 36 -bS | \
	samtools sort -@ 36 -o $datadir/rna_align/nivar_fb15_cDNA.sorted.bam -T $datadir/rna_align/reads.tmp -
    samtools index $datadir/rna_align/nivar_fb15_cDNA.sorted.bam
fi


if [ $1 == get_rna_regions ] ; then
    datadir=/kyber/Data/seqlab/sp_2019/nivar_r9

    bedtools genomecov -bg -ibam $datadir/rna_align/nivar_fb15_cDNA.sorted.bam | awk '$4 > 4' | bedtools merge -d 2 > $datadir/rna_align/nivar_fb15_cDNA.bed
    bedtools genomecov -bg -ibam $datadir/rna_align/nivar_fb15_dRNA.sorted.bam | awk '$4 > 4' | bedtools merge -d 2 > $datadir/rna_align/nivar_fb15_dRNA.bed
fi

if [ $1 == find_intersects ] ; then
    datadir=/kyber/Data/seqlab/sp_2019/nivar_r9
    anndir=$datadir/annotation
    
    bedtools intersect -wa -a $anndir/augustus/cani_by_caal.gff -b $datadir/rna_align/nivar_fb15_dRNA.bed > $datadir/rna_align/nivar_fb15_dRNA_genes.bed
    bedtools intersect -wa -a $anndir/augustus/cani_by_caal.gff -b $datadir/rna_align/nivar_fb15_cDNA.bed > $datadir/rna_align/nivar_fb15_cDNA_genes.bed
fi

 
     
