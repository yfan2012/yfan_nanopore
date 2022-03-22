#!/bin/bash

projdir=/mithril/Data/Nanopore/projects/methbin
datadir=$projdir/barnyard/strains

if [ $1 == pick_dcm_loci ] ; then
    mkdir -p $datadir/modcheck

    emegafile=$datadir/megalodon/ecoli_plas/per_read_modified_base_calls.txt
    smegafile=$datadir/megalodon/staph_plas/per_read_modified_base_calls.txt
    
    ##ecoli chr
    awk '($4=="4308216" || $4=="4381955" || $4=="484397" || $4=="166504" || $4=="2201921" || $4=="4308218" || $4=="4381957" || $4=="484399" || $4=="166506" || $4=="2201923") && $2=="CP017100.1" {print $0}' $emegafile > $datadir/modcheck/ecoli_dcm.txt &
    
    ##ecoli plas
    awk '($4=="2264" || $4=="7289" || $4=="8577" || $4=="2266" || $4=="7291" || $4=="8579") && $2=="PRW62" {print $0}' $emegafile > $datadir/modcheck/ecoliplas_dcm.txt &

    ##staph chr
    awk '($4=="1289381" || $4=="772611" || $4=="1766568" || $4=="548492" || $4=="292706" || $4=="1289383" || $4=="772613" || $4=="1766570" || $4=="548494" || $4=="292708") && $2=="NC_007795.1"  {print $0}' $smegafile > $datadir/modcheck/staph_dcm.txt &

    ##staph plas
    awk '($4=="2264" || $4=="7289" || $4=="8577" || $4=="2266" || $4=="7291" || $4=="8579") && $2=="PRW62" {print $0}' $smegafile > $datadir/modcheck/staphplas_dcm.txt
fi


if [ $1 == pick_dam_loci ] ; then
    emegafile=$datadir/megalodon/ecoli_plas/per_read_modified_base_calls.txt
    smegafile=$datadir/megalodon/staph_plas/per_read_modified_base_calls.txt

    
    ##ecoli chr
    awk '($4=="563109" || $4=="895635" || $4=="1247771" || $4=="2907090" || $4=="328450" || $4=="563110" || $4=="895636" || $4=="1247772" || $4=="2907091" || $4=="328451") && $2=="CP017100.1" {print $0}' $emegafile > $datadir/modcheck/ecoli_dam.txt &
    
    ##ecoli plas
    awk '($4=="5411" || $4=="1706" || $4=="3746" || $4=="7432" || $4=="8341" || $4=="5412" || $4=="1707" || $4=="3747" || $4=="7433" || $4=="8342") && $2=="PRW62" {print $0}' $emegafile > $datadir/modcheck/ecoliplas_dam.txt &

    ##staph chr
    awk '($4=="922259" || $4=="1922606" || $4=="832844" || $4=="2726130" || $4=="2497086" || $4=="922260" || $4=="1922607" || $4=="832845" || $4=="2726131" || $4=="2497087") && $2=="NC_007795.1"  {print $0}' $smegafile > $datadir/modcheck/staph_dam.txt &

    ##staph plas
    awk '($4=="5411" || $4=="1706" || $4=="3746" || $4=="7432" || $4=="8341" || $4=="5412" || $4=="1707" || $4=="3747" || $4=="7433" || $4=="8342") && $2=="PRW62" {print $0}' $smegafile > $datadir/modcheck/staphplas_dam.txt &
    
fi
