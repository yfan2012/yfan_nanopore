#!/bin/bash


seqdir=~/Dropbox/Lab/carbapenem_r21/annotations
i=$seqdir/171002_bucc.spades.contigs.fasta
outdir=~/Dropbox/Lab/carbapenem_r21/annotations


file=`echo $i | cut -d / -f 8`
samp=`echo $file | cut -d . -f 1,2`


abriout=$outdir/abricate
abricate $seqdir/$file > $abriout/$samp.raw.tsv

prokout=$outdir/prokka
prokka --outdir $prokout --genus Klebsiella --usegenus --prefix ${samp}_prokka --cpus 12 --force $i &> $prokout/log/$samp.prokka_log.txt

rgiout=$outdir/rgicard
rgi -t contig -i $i -n 8 -o $rgiout/$samp.rgicard


