#!/bin/bash

prefixes='5_millions_reads_psample_APCflfl_Exp1 1_7_million_psampleAPCflfl_Exp1'

if [ $1 == gatherfq ] ; then
    for prefix in $prefixes ;
    do
	bash abundance_analysis.sh gatherfq $prefix
    done
fi

if [ $1 == ontkraken ] ; then
    for prefix in $prefixes ;
    do
	bash abundance_analysis.sh ontkraken $prefix
    done
fi

if [ $1 == stdkraken ] ; then
    for prefix in $prefixes ;
    do
	bash abundance_analysis.sh stdkraken $prefix
    done
fi

if [ $1 == ontbracken ] ; then
    for prefix in $prefixes ;
    do
	echo %prefix
	bash abundance_analysis.sh ontbracken $prefix
    done
fi

if [ $1 == stdbracken ] ; then
    for prefix in $prefixes ;
    do
	bash abundance_analysis.sh stdbracken $prefix
    done
fi

