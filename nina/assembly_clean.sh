#!/bin/bash

datadir=~/Dropbox/yfan/nina_fungus/assemblies

if [ $1 == mum_assemblers ] ; then
    ##mummer assemblies against each other. Try to ditch wtdbg2 junk contigs in some way
    for i in st31 st90853 ;
    do
	
    
