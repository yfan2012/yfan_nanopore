#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/dunlop/180714_dunlop_3ecoli
srcdir=~/Code/utils/marcc

##mkdir -p $datadir/batch_logs

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_dunlop $srcdir/untar.scr $datadir/180714_dunlop_3ecoli.tar.gz $datadir
fi

if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_logs
    mkdir -p $datadir/call_done
    sbatch --array=0-693 --job-name=call_dunlop --output=$datadir/call_logs/dunlop_ecoli3.%A_%a.out $srcdir/bc_call.scr $datadir
fi


if [ $1 == fastq ] ; then
    for i in 1 2 3 ;
    do
	mkdir -p $datadir/ecoli${i}
	cat $datadir/called/*/workspace/pass/barcode0${i}/*fastq > $datadir/ecoli${i}/fastqs/ecoli${i}.fastq
    done
fi
    
if [ $1 == downsamp ] ; then
    for i in 1 2 3 ;
    do
	head -800000 $datadir/ecoli${i}/fastqs/ecoli${i}.fastq > $datadir/ecoli${i}/fastqs/ecoli${i}_sub200k.fastq
    done
fi

if [ $1 == assemble ] ; then
    for i in 1 2 3 ;
    do
	mkdir $datadir/ecoli${i}/assembly
	bash assemble_bacteria.sh $datadir/ecoli${i}/fastqs/ecoli${i}_sub200k.fastq $datadir/ecoli${i}/assembly
    done
fi

if [ $1 == illumina_dilith ] ; then
    ##copy illumina data to dilithium
    mkdir -p /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli
    for i in `find ~/BaseSpace/Projects/180722_dunlop_3ecoli/ -name *fastq.gz` ;
    do
	cp $i /dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/
    done
fi

if [ $1 == illumina_marcc ] ; then
    mkdir -p $datadir/illumina
    scp -r smaug:/dilithium/Data/NGS/Raw/180722_dunlop_3ecoli/* $datadir/illumina/
fi

if [ $1 == pilon ] ; then
    for i in 1 2 3 ;
    do
	cp $datadir/illumina/ecoli${i}*.fastq.gz $datadir/ecoli${i}/pilon/
	sbatch --output=$datadir/batch_logs/ecoli${i}_pilon.out --job-name=ecoli${i} pilon.scr $datadir/ecoli${i}	
    done
fi


if [ $1 == copyassembly ] ; then
    scp -r $datadir/ecoli* smaug:/dilithium/Data/Nanopore/projects/dunlop_ecoli/
fi


if [ $1 == parsnp ] ; then
    parsdir=~/Dropbox/yfan/dunlop/parsnp
    assembledir=~/Dropbox/yfan/dunlop/assemblies
    
    mkdir -p $parsdir/raw_parsnp
    mkdir -p $parsdir/canu17_parsnp
    mkdir -p $parsdir/pilon_parsnp
    mkdir -p $parsdir/pilon17_parsnp

    parsnp -r $assembledir/bw25311.fasta -d $assembledir/raw -p 12 -o $parsdir/raw_parsnp -c
    harvesttools -i $parsdir/raw_parsnp/parsnp.ggr -V $parsdir/raw_parsnp/raw.vcf
    parsnp -r $assembledir/bw25311.fasta -d $assembledir/canu17 -p 12 -o $parsdir/canu17_parsnp -c
    harvesttools -i $parsdir/canu17_parsnp/parsnp.ggr -V $parsdir/canu17_parsnp/canu17.vcf
    parsnp -r $assembledir/bw25311.fasta -d $assembledir/pilon -p 12 -o $parsdir/pilon_parsnp -c
    harvesttools -i $parsdir/pilon_parsnp/parsnp.ggr -V $parsdir/pilon_parsnp/pilon.vcf
    parsnp -r $assembledir/bw25311.fasta -d $assembledir/pilon17 -p 12 -o $parsdir/pilon17_parsnp -c
    harvesttools -i $parsdir/pilon17_parsnp/parsnp.ggr -V $parsdir/pilon17_parsnp/pilon17.vcf
fi

if [ $1 == spades ] ; then
    ml python/2.7
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	fastqdir=~/work/180714_dunlop_3ecoli/illumina
	sampdir=$datadir/${i}_spades
	mkdir $sampdir
	mkdir $fastqdir/trimmed
	
	java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $fastqdir/$i*R1* $fastqdir/$i*R2* $fastqdir/trimmed/${i}_forward_paired.fq.gz $fastqdir/trimmed/${i}_forward_unpaired.fq.gz $fastqdir/trimmed/${i}_reverse_paired.fq.gz $fastqdir/trimmed/${i}_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	
	spades.py -1 $fastqdir/trimmed/${i}_forward_paired.fq.gz -2 $fastqdir/trimmed/${i}_reverse_paired.fq.gz -t 36 -m 300 -o $datadir/${i}_spades
    done
fi



if [ $1 == reorg ] ; then
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	mv -p $datadir/$i
	mv $datadir/${i}_* $datadir/$i
	mv $datadir/$i/*assembly $datadir/$i/assembly
	mv $datadir/$i/*assembly17 $datadir/$i/assembly17
	mv $datadir/$i/*pilon $datadir/$i/pilon
	mv $datadir/$i/*pilon17 $datadir/$i/pilon17
    done
fi
	
if [ $1 == align_pilon17diff ] ; then
    ##look at persistent CCWGG snp in ecoli 3 by alignment
    ml samtools
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	fastqdir=~/work/180714_dunlop_3ecoli/illumina
	bamdir=$datadir/$i/align
	pilonill=$datadir/$i/align/pilonill
	
	mkdir -p $bamdir
	mkdir -p $pilonill
	mkdir -p $pilonill/btidx
	
	cp $datadir/$i/pilon17/$i.pilon17.10.fasta $pilonill/btidx/
	bowtie2-build -q $pilonill/btidx/$i.pilon17.10.fasta $pilonill/btidx/ecoli3.pilon17.10
	bowtie2 -p 24 -x $pilonill/btidx/ecoli3.pilon17.10 -1 $datadir/illumina/$i*R1_001.fastq.gz -2 $datadir/illumina/$i*R2_001.fastq.gz | samtools view -bS - | samtools sort -o $pilonill/$i.sorted.bam
	samtools index $pilonill/$i.sorted.bam
    done
fi

	
if [ $1 == cppilonill ] ; then
   for i in ecoli1 ecoli2 ecoli3 ;
   do
       pilonill=$datadir/$i/align/pilonill
       scp -r $pilonill/$i.sorted.bam* smaug:~/Dropbox/yfan/dunlop/align/pilonill/
   done
fi
       
if [ $1 == align_pilondiff ] ; then
    ##look at persistent CCWGG snp in ecoli 3 by alignment
    ml samtools
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	fastqdir=~/work/180714_dunlop_3ecoli/illumina
	bamdir=$datadir/$i/align
	pilonill=$datadir/$i/align/pilonill
	
	mkdir -p $bamdir
	mkdir -p $pilonill
	mkdir -p $pilonill/btidx
	
	cp $datadir/$i/pilon/$i.pilon.10.fasta $pilonill/btidx/
	bowtie2-build -q $pilonill/btidx/$i.pilon.10.fasta $pilonill/btidx/ecoli3.pilon.10
	bowtie2 -p 24 -x $pilonill/btidx/ecoli3.pilon.10 -1 $datadir/illumina/$i*R1_001.fastq.gz -2 $datadir/illumina/$i*R2_001.fastq.gz | samtools view -bS - | samtools sort -o $pilonill/$i.sorted.bam
	samtools index $pilonill/$i.sorted.bam
    done
fi


if [ $1 == kleb_grant ] ; then
    ##mummer klpn_70 to show insertion element
    datadir=~/Dropbox/Timplab_Data/cpowgs/mummer/KLPN_70
    asmdir=~/Dropbox/Timplab_Data/cpowgs/assemblies/pilon
    
    nucmer -p $datadir/KLPN_70 $asmdir/KLPN_139.pilon.10.fasta $asmdir/KLPN_70.pilon.10.fasta 
    mummerplot --png -p $datadir/KLPN_70.layout $datadir/KLPN_70.delta -R $asmdir/KLPN_139.pilon.10.fasta -Q $asmdir/KLPN_70.pilon.10.fasta
    mummerplot --png -p $datadir/KLPN_70 $datadir/KLPN_70.delta
fi


if [ $1 == kleb_mgrb ] ; then
    srcdir=~/Code/carbapenem_r21/mutations
        
    datadir=/atium/Data/Nanopore/cpowgs
    contigs=$datadir/assemblies/spades/KLPN_70.spades.fasta
    genedir=$datadir/References/ref_gene_seqs
    genes=$genedir/klpn_genes.fasta
    
    blastdbdir=$datadir/blastdb
    qual=spades
    name=KLPN_70
    
    mkdir -p $blastdbdir/$qual
    mkdir -p $datadir/blast_hits/$qual
    mkdir -p $datadir/geneseqs/$qual
    
    makeblastdb -in $contigs -out $blastdbdir/$qual/$name.db -dbtype nucl
    blastn -query $genedir/klpn_genes.fasta -db $blastdbdir/$qual/$name.db -outfmt 7 -out $datadir/blast_hits/$qual/$name.tsv
    python $srcdir/make_gene_fa.py -t $datadir/blast_hits/$qual/$name.tsv -a $contigs -g $genes -o $datadir/geneseqs/$qual > $datadir/geneseqs/$qual/$name.short.tsv
    Rscript $srcdir/plot_shortalign.R -s $datadir/geneseqs/$qual/$name.short.tsv -o ~/Dropbox/yfan/carbapenem_r21/mutations/gene_cov/KLPN_70_spades.pdf
    Rscript $srcdir/plot_shortalign.R -s $datadir/geneseqs/polished/$name.short.tsv -o ~/Dropbox/yfan/carbapenem_r21/mutations/gene_cov/KLPN_70.pdf
fi
   
