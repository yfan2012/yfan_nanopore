#!/bin/bash

datadir=~/work/r10_zymo

if [ $1 == assembly ] ; then
    mkdir $datadir/canu_r10
    mkdir $datadir/canu_r9
    
    canu \
	-p zymo_r10 -d $datadir/canu_r10/ \
	-gridOptions="--time=22:00:00 --partition=parallel" \
	genomeSize=15m \
	corMinCoverage=0 \
	corOutCoverage=all \
	corMhapSensitivity=high \
	corMaxEvidenceCoverageLocal=10 \
	corMaxEvidenceCoverageGlobal=10 \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/fastqs/zymo_r10.fq.gz
    canu \
	-p zymo_r9 -d $datadir/canu_r9/ \
	-gridOptions="--time=22:00:00 --partition=parallel" \
	genomeSize=15m \
	corMinCoverage=0 \
	corOutCoverage=all \
	corMhapSensitivity=high \
	corMaxEvidenceCoverageLocal=10 \
	corMaxEvidenceCoverageGlobal=10 \
	stopOnReadQuality=false \
	-nanopore-raw $datadir/fastqs/zymo_r9.fq.gz
fi
