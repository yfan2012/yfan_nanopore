#!/bin/bash

datadir=/atium/Data/Nanopore/cpowgs
seqsdir=$datadir/References/gene_seqs
ref=/mithril/Data/NGS/Reference/kpneumo/NC_016845.fasta

##make the database out of the reference genome
makeblastdb -in $ref -out $datadir/blastdb/refdb -dbtype nucl

##combine gene sequences into one fasta
for i in $seqsdir/*fasta;
do
    name=`echo $i | cut -d '/' -f 8 | cut -d '.' -f 1 `
    fasthead='>'${name}
    sed -i "1s/.*/$fasthead/" $i
done

cat $seqsdir/*fasta > $seqsdir/all_genes.fa


##Search for genes in reference
blastn -query $seqsdir/all_genes.fa -db $datadir/blastdb/refdb -outfmt 7 -out $datadir/genomes/var_profiles/ref_genes.tsv
blastn -query $seqsdir/all_genes.fa -subject $ref -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines
blastn -query $seqsdir/all_genes.fa -subject $ref -dust no -parse_deflines -out $datadir/genomes/var_profiles/ref_genes.vis


##Pull out reference genes based on above
##python make_ref_genes.py


##combine
refgenesdir=$datadir/References/ref_gene_seqs
cat $refgenesdir/*fasta > $refgenesdir/all_genes.fa


##Make database out of all assemblies with raw and pilon assemblies

i=170816_BUCC
pilon=$datadir/$i/pilon/$i.fasta
makeblastdb -in $pilon -out $datadir/blastdb/$i.pilondb -dbtype nucl
##blastn -query $refgenesdir/all_genes.fa -db $datadir/blastdb/$i.pilondb -outfmt 7 -out $datadir/genomes/var_profiles/$i.pilon_genes.tsv
##blastn -query $refgenesdir/all_genes.fa -subject $pilon -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$i.pilon_genes.btop
##blastn -query $refgenesdir/all_genes.fa -subject $pilon -dust no -parse_deflines -out $datadir/genomes/var_profiles/$i.pilon_genes.vis

blastn -query $seqsdir/all_genes.fa -db $datadir/blastdb/$i.pilondb -outfmt 7 -out $datadir/genomes/var_profiles/$i.pilon_genes_orig.tsv
blastn -query $seqsdir/all_genes.fa -subject $pilon -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$i.pilon_genes_orig.btop
blastn -query $seqsdir/all_genes.fa -subject $pilon -dust no -parse_deflines -out $datadir/genomes/var_profiles/$i.pilon_genes_orig.vis


canu=$datadir/$i/canu_assembly/$i.contigs.fasta
makeblastdb -in $canu -out $datadir/blastdb/$i.rawdb -dbtype nucl
##blastn -query $refgenesdir/all_genes.fa -db $datadir/blastdb/$i.rawdb -outfmt 7 -out $datadir/genomes/var_profiles/$i.raw_genes.tsv
##blastn -query $refgenesdir/all_genes.fa -subject $canu -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$i.raw_genes.btop
##blastn -query $refgenesdir/all_genes.fa -subject $canu -dust no -parse_deflines -out $datadir/genomes/var_profiles/$i.raw_genes.vis

blastn -query $seqsdir/all_genes.fa -db $datadir/blastdb/$i.rawdb -outfmt 7 -out $datadir/genomes/var_profiles/$i.raw_genes_orig.tsv
blastn -query $seqsdir/all_genes.fa -subject $canu -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$i.raw_genes_orig.btop
blastn -query $seqsdir/all_genes.fa -subject $canu -dust no -parse_deflines -out $datadir/genomes/var_profiles/$i.raw_genes_orig.vis


illumina=$datadir/$i/171002_bucc.spades.contigs.fasta
makeblastdb -in $illumina -out $datadir/blastdb/$i.illuminadb -dbtype nucl
##blastn -query $refgenesdir/all_genes.fa -db $datadir/blastdb/$i.illuminadb -outfmt 7 -out $datadir/genomes/var_profiles/$i.illumina_genes.tsv
##blastn -query $refgenesdir/all_genes.fa -subject $illumina -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$i.illumina_genes.btop
##blastn -query $refgenesdir/all_genes.fa -subject $illumina -dust no -parse_deflines -out $datadir/genomes/var_profiles/$i.illumina_genes.vis

blastn -query $seqsdir/all_genes.fa -db $datadir/blastdb/$i.illuminadb -outfmt 7 -out $datadir/genomes/var_profiles/$i.illumina_genes_orig.tsv
blastn -query $seqsdir/all_genes.fa -subject $illumina -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $datadir/genomes/var_profiles/$i.illumina_genes_orig.btop
blastn -query $seqsdir/all_genes.fa -subject $illumina -dust no -parse_deflines -out $datadir/genomes/var_profiles/$i.illumina_genes_orig.vis


    
