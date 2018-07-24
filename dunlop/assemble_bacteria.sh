#!/bin/bash

prefix=`basename $1 .fastq`

canu \
    -p $prefix -d $2 \
    -gridOptions="--mem=8g --time=22:00:00 --account=mschatz1" \
    genomeSize=4.5m \
    -nanopore-raw $1
