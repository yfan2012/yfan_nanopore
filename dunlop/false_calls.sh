#!bin/bash

datadir=~/Dropbox/yfan/dunlop
asmdir=$datadir/assemblies
alndir=$datadir/align

if [ $1 == mummer ] ; then
    nucdir=$datadir/mummer
    mkdir -p $nucdir
    for i in ecoli1 ecoli2 ecoli3 ;
    do
	dnadiff -p $nucdir/$i.pilon $asmdir/bw25311.fasta $asmdir/pilon/$i.pilon.10.fasta
	dnadiff -p $nucdir/$i.pilon17 $asmdir/bw25311.fasta $asmdir/pilon17/$i.pilon17.10.fasta
    done
fi
