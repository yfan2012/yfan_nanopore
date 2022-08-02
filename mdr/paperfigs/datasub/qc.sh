#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin
zymodir=$datadir/zymo
barndir=$datadir/barnyard/strains
mdrdir=$datadir/mdr

dbxdir=~/gdrive/mdr/paperfigs/qc

if [ $1 == yield ] ; then
    zymofq=$zymodir/fastq/20190809_zymo_control/20190809_zymo_control.fq.gz
    staphfq=$barndir/fastqs/220131_mdr_barnyard_st3294.fastq.gz
    ecolifq=$barndir/fastqs/220131_mdr_barnyard_st3689.fastq.gz
    mdrfq=$mdrdir/fastqs/200708_mdr_stool16native_nohuman.fq.gz

    phase1=$mdrdir/dataclean/181127_hiC_stool_nohuman/181127_hiC_stool_phase.humanfilt.fastq_1.fastq.gz
    phase2=$mdrdir/dataclean/181127_hiC_stool_nohuman/181127_hiC_stool_phase.humanfilt.fastq_2.fastq.gz
    shot1=$mdrdir/dataclean/181127_hiC_stool_nohuman/181127_hiC_stool_shotgun.humanfilt.fastq_1.fastq.gz
    shot2=$mdrdir/dataclean/181127_hiC_stool_nohuman/181127_hiC_stool_shotgun.humanfilt.fastq_2.fastq.gz
    
    touch $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $zymofq >> $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $staphfq >> $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $ecolifq >> $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $mdrfq >> $dbxdir/yields.csv
    
    bash ~/Code/utils/qc/basic_run_assess.sh $phase1 >> $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $phase2 >> $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $shot1 >> $dbxdir/yields.csv
    bash ~/Code/utils/qc/basic_run_assess.sh $shot2 >> $dbxdir/yields.csv

fi
