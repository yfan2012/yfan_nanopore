#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/nivar
datadir=/uru/Data/Nanopore/projects/nivar/paperfigs
fq=/uru/Data/Nanopore/projects/nivar/r9/r9_3kb.fq
dbxdir=~/Dropbox/yfan/nivar/paperfigs

ref=$rawdir/reference/candida_nivariensis.fa
gla=$rawdir/reference/medusa_fungi/candida_glabrata.fa

fin=$datadir/assembly_final/nivar.final.fasta
gff=$datadir/annotation_final/nivar.final.gff

xugla=$datadir/patho/glabrata_xu.fa

if [ $1 == amino ] ; then
    gffread \
	-y $datadir/patho/long_genes.fa \
	-g $fin \
	$datadir/patho/long_genes.gff
fi

if [ $1 == subtelo ] ; then
    ##blast ref and asm for list of subtelo genes from xu glabrata study
    makeblastdb \
	-in $fin \
	-out $datadir/patho/nivar.final \
	-dbtype nucl

    blastn \
	-query $datadir/patho/gpicwp.fa \
	-db $datadir/patho/nivar.final \
	-outfmt 7 \
	-out $datadir/patho/nivar.final.gpicwp_hits.tsv

    makeblastdb \
	-in $ref \
	-out $datadir/patho/candida_nivariensis \
	-dbtype nucl

    blastn \
	-query $datadir/patho/gpicwp.fa \
	-db $datadir/patho/candida_nivariensis \
	-outfmt 7 \
	-out $datadir/patho/candida_nivariensis.tsv
fi

if [ $1 == adhesion ] ; then
    ##blast ref and asm for list of putative adhesion prots from predGPI
    makeblastdb \
	-in $fin \
	-out $datadir/patho/nivar.final \
	-dbtype nucl

    blastn \
	-query $datadir/patho/gpigenes.fasta \
	-db $datadir/patho/nivar.final \
	-outfmt 7 \
	-out $datadir/patho/nivar.final.predgpi_hits.tsv

    makeblastdb \
	-in $ref \
	-out $datadir/patho/candida_nivariensis \
	-dbtype nucl

    blastn \
	-query $datadir/patho/gpigenes.fasta \
	-db $datadir/patho/candida_nivariensis \
	-outfmt 7 \
	-out $datadir/patho/candida_nivariensis.predgpi_hits.tsv
fi
