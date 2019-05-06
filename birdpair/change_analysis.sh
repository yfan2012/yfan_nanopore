#!/bin/bash

datadir=~/Dropbox/timplab_data/birdpair
asmdir=$datadir/assemblies/canu

if [ $1 == cleanfiles ] ; then
    sed -i -e 's/_pilon//g' $asmdir/bird_canu.pilon.10*
    sed -i -e 's/>/>pilon_/g' $asmdir/bird_canu.pilon.10*
    sed -i -e 's/_pilon//g' $asmdir/patient_canu.pilon.9*
    sed -i -e 's/>/>pilon_/g' $asmdir/patient_canu.pilon.9*
fi

if [ $1 == mummer ] ; then
    ##do all the mummer things
    mkdir -p $datadir/mummer
    mkdir -p ~/tmp

    cp $asmdir/bird_canu.pilon.10.fasta ~/tmp/
    cp $asmdir/patient_canu.pilon.9.fasta ~/tmp/
    
    nucmer -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    mummerplot --fat --png -p ~/tmp/birdpair.layout ~/tmp/birdpair.delta -R ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta
    mummerplot --fat --png -p ~/tmp/birdpair ~/tmp/birdpair.delta

    dnadiff -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    cp ~/tmp/birdpair* $datadir/mummer/
fi

if [ $1 == mummerfilt ] ; then
    ##do all the mummer things
    mkdir -p $datadir/mummer_filter
    mkdir -p ~/tmp

    cp $asmdir/bird_canu.pilon.10.fasta ~/tmp/
    cp $asmdir/patient_canu.pilon.9.fasta ~/tmp/
    
    nucmer -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    mummerplot --filter --fat --png -p ~/tmp/birdpair.layout ~/tmp/birdpair.delta -R ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta
    mummerplot --filter --fat --png -p ~/tmp/birdpair ~/tmp/birdpair.delta

    dnadiff -p ~/tmp/birdpair ~/tmp/bird_canu.pilon.10.fasta ~/tmp/patient_canu.pilon.9.fasta

    cp ~/tmp/birdpair* $datadir/mummer_filter/
fi
