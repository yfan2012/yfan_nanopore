#!/bin/bash

datadir=/dilithium/Data/Nanopore/sindbis

for i in antibody infected mock Antibody_2dpi Antibody_3dpi Sindbis_2dpi Sindbis_3dpi ;
do
    if [ $1 == stringtie ] ; then
	mkdir -p $datadir/$i/stringtie
	stringtie -p 36 -L -G $datadir/rattus_norvegicus.gff -o $datadir/$i/stringtie/$i.rat.gtf $datadir/$i/align/$i.rat.splicealn.sorted.bam
    fi
done
       
