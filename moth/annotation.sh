#!/bin/bash

datadir=/uru/Data/Nanopore/projects/moth
asm=$datadir/correction2_moth_scaffolds.fasta

if [ $1 == augustus ] ; then
    ##try with all bugs
    ##adorsata=giant honeybee
    ##aedes=mosquito
    ##bombus=bumblebee
    ##camponotus_floridanus=carpenter ant
    ##culex=mosquito
    ##nasonia=wasp
    for i in adorsata aedes bombus_impatiens1 bombus_terrestris2 camponotus_floridanus culex fly nasonia rhodnius ;
    do
	augustus --species=$i $asm &> $datadir/moth_annotation.$i.gff &
    done
fi


