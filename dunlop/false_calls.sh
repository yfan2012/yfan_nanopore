#!bin/bash

datadir=~/Dropbox/Timplab_Data/dunlop
asmdir=$datadir/assemblies
alndir=$datadir/align

if [ $1 == mummer ] ; then
    nucdir=$datadir/mummer
    mkdir -p $nucdir
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	nucmer -p $nucdir/$i $asmdir/bw25311.fasta $asmdir/ecoli3/pilon17/$i.pilon17.10.fasta
	mummerplot --png -p $nucdir/$i.layout $nucdir/$i.delta -R $asmdir/bw25311.fasta -Q $asmdir/ecoli3/pilon17/$i.pilon17.10.fasta
	mummerplot --png -p $nucdir/$i $nucdir/$i.delta

	##dnadiff -p $nucdir/$i.pilon $asmdir/bw25311.fasta $asmdir/pilon/$i.pilon.10.fasta
	##dnadiff -p $nucdir/$i.pilon17 $asmdir/bw25311.fasta $asmdir/pilon17/$i.pilon17.10.fasta
    done
fi
