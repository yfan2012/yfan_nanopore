#!/bin/bash

datadir=/pym/Data/Nanopore/projects/prolificans

if [ $1 == mummer ] ; then
    for i in st31 st90853 st5317 ;
    do
	mkdir -p $datadir/$i/mummer

	
