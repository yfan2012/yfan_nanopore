#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr
dbxdir=~/Dropbox/timplab_data/mdr

if [ $1 == abricate ] ; then
    ##detect genes in all metagenomic asms
    mkdir -p $datadir/MDRstool_19/amr

    for db in resfinder argannot card ecoh ncbi plasmidfinder vfdb megares ecoli_vf ;
    do

	for i in 100m 10g 150m 500m ;
	do
	    asm=$datadir/MDRstool_19/metaflye/$i/MDRstool_19_$i.assembly.fasta
	    prefix=`basename $asm .assembly.fasta`
	    abricate --quiet --threads 36 --db $db --nopath $asm > $datadir/MDRstool_19/amr/$prefix.$db.tsv
	done

	illasm=$datadir/illumina/metaspades/181127_hiC_stool_shotgun.contigs.fasta
	abricate --quiet --threads 36 --db $db --nopath $illasm > $datadir/MDRstool_19/amr/MDRstool_19_metaspades.$db.tsv
    done
fi


if [ $1 == stool16_abricate ] ; then
    ##detect genes in mdr16, since it's now unclear which one has illumina data
    mkdir -p $datadir/MDRstool_16/amr

    for db in resfinder argannot card ecoh ncbi plasmidfinder vfdb megares ecoli_vf ;
    do

	for i in 100m 10g 150m 500m ;
	do
	    asm=$datadir/MDRstool_16/metaflye/$i/MDRstool_16_$i.assembly.fasta
	    prefix=`basename $asm .assembly.fasta`
	    abricate --quiet --threads 36 --db $db --nopath $asm > $datadir/MDRstool_16/amr/$prefix.$db.tsv
	done
    done
fi


if [ $1 == cat_abricate ] ; then
    ##merge amr files

    ##for i in MDRstool_19 MDRstool_16 illumina ;
    for i in MDRstool_19 MDRstool_16 ;
    do
	for samp in $datadir/$i/amr/*card*.tsv ;
	do
	    prefix=` basename $samp .card.tsv `
	    sed -s 1d $datadir/$i/amr/${prefix}*.tsv > $datadir/$i/amr/${prefix}.all.tsv
	done
    done
fi

if [ $1 == merge_all ] ; then
    ##merge samples

    stool19=$datadir/MDRstool_19/amr/*all.tsv
    stool16=$datadir/MDRstool_16/amr/*all.tsv
    illumina=$datadir/illumina/amr/*all.tsv

    header=`head -n 1 $datadir/MDRstool_19/amr/MDRstool_19_100m.card.tsv`

    sed -i -e "1i $header" $stool19
    sed -i -e "1i $header" $stool16
    sed -i -e "1i $header" $illumina
    
    abricate --summary $stool19 $stool16 $illumina > $datadir/all/abricate_all.tsv
fi


if [ $1 == combine_db_fastas ] ; then
    ##make some phylo trees of all genes that abricate looks at, and maybe this will help group them?
    cat ~/software/abricate/db/*/sequences > $datadir/all/abricate_alldb.fasta
fi

if [ $1 == try_clustalw ] ; then

    clustalw -i $datadir/all/abricate_alldb.fasta \
	     -o $datadir/all/clustalw.fasta \
	     --distmat-out $datadir/all/clustalw_dists.txt \
	     --percent-id \
	     --threads 54 \
	     --force \
	     --full \
	     -t DNA \
	     -v
fi

if [ $1 == convertcsv ] ; then
    sed -e 's/\s\+/,/g' $datadir/amr/clustalw_dists.txt > $datadir/amr/clustalw_dists.csv
fi


if [ $1 == new_stool16_abricate ] ; then
    ##detect genes in mdr16, since it's now unclear which one has illumina data
    mkdir -p $datadir/MDRstool_16/amr

    for db in resfinder argannot card ecoh ncbi plasmidfinder vfdb megares ecoli_vf ;
    do

	for i in native_100m pcr_100m ;
	do
	    asm=$datadir/MDRstool_16/metaflye/$i/MDRstool_16_$i.assembly.fasta
	    prefix=`basename $asm .assembly.fasta`
	    abricate --quiet --threads 36 --db $db --nopath $asm > $datadir/MDRstool_16/amr/$prefix.$db.tsv
	done
    done
fi

if [ $1 == new_merge_all ] ; then
    ##merge samples

    native=$datadir/MDRstool_16/amr/*native*all.tsv
    pcr=$datadir/MDRstool_16/amr/*pcr*all.tsv
    illumina=$datadir/illumina/amr/*all.tsv

    header=`head -n 1 $datadir/MDRstool_16/amr/MDRstool_16_native_100m.card.tsv`

    sed -i -e "1i $header" $native
    sed -i -e "1i $header" $pcr
    
    abricate --summary $pcr $native $illumina > $datadir/all/abricate_all.tsv
fi
