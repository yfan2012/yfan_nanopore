#!/bin/bash

##sequencing/basecalling done on the grid

datadir=/dilithium/Data/Nanopore/sindbis
ref=

if [ $1 == align ] ; then
    
    minimap2 -ax splice -uf -k14
