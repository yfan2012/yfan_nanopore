#!/bin/bash

##make blast db from all penA genes in bcc
##blast hometown burkle against database
##win

seqsdir=/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_penA_genes
outdir=~/Dropbox/Lab/carbapenem_r21/annotations/bucc/ident

##hometown burkle as asked for by p-dawg
htburk=$seqsdir/penA_bcepacia_170816_BUCC_pilon.fasta
makeblastdb -in $htburk -out $seqsdir/penA_bcepacia_170816_BUCC_pilondb -dbtype nucl

blastn -query $seqsdir/bcc_penA_all.fa -db $seqsdir/penA_bcepacia_170816_BUCC_pilondb -outfmt 7 -out $outdir/htburk_identity.tsv
blastn -query $seqsdir/bcc_penA_all.fa -subject $htburk -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $outdir/htburk_identity.btop
blastn -query $seqsdir/bcc_penA_all.fa -subject $htburk -dust no -parse_deflines -out $outdir/htburk_identity.vis


##contaminans as requested by t-dawg
contam1=$seqsdir/penA_bcontaminans1.fasta
makeblastdb -in $contam1 -out $seqsdir/penA_bcontaminans1db -dbtype nucl

blastn -query $seqsdir/bcc_penA_all.fa -db $seqsdir/penA_bcontaminans1db -outfmt 7 -out $outdir/contam1_identity.tsv
blastn -query $seqsdir/bcc_penA_all.fa -subject $contam1 -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $outdir/contam1_identity.btop
blastn -query $seqsdir/bcc_penA_all.fa -subject $contam1 -dust no -parse_deflines -out $outdir/contam1_identity.vis


contam2=$seqsdir/penA_bcontaminans2.fasta
makeblastdb -in $contam2 -out $seqsdir/penA_bcontaminans2db -dbtype nucl

blastn -query $seqsdir/bcc_penA_all.fa -db $seqsdir/penA_bcontaminans2db -outfmt 7 -out $outdir/contam2_identity.tsv
blastn -query $seqsdir/bcc_penA_all.fa -subject $contam2 -dust no -outfmt "6 qseqid sseqid btop" -parse_deflines -out $outdir/contam2_identity.btop
blastn -query $seqsdir/bcc_penA_all.fa -subject $contam2 -dust no -parse_deflines -out $outdir/contam2_identity.vis
