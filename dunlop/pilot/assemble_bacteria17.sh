#!/bin/bash

prefix=` basename $1 _sub200k.fastq `
echo $prefix


##Assemble if it's not already done

~/software/canu-1.7/Linux-amd64/bin/canu \
    -p $prefix -d $2 \
    -gridOptions="--mem=8g --time=22:00:00 --account=mschatz1" \
    genomeSize=4.5m \
    -nanopore-raw $1

