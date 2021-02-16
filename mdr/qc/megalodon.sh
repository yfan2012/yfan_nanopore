#!/bin/bash

##wrap for easy running of megalodon with how I've set up the neb data

samp=$1 ##name of neb samp
##motif=$2 ##motif


ref=/mithril/Data/Nanopore/projects/methbin/reference/allsamps.fa

if [ $2 == rerio_allmod ] ; then
	mkdir -p ~/data/mdr/qc/megalodon/$samp/${samp}

	time megalodon \
	    ~/data/mdr/qc/multiraw_sub/$samp \
	    --overwrite \
	    --guppy-server-path "/usr/bin/guppy_basecall_server" \
	    --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	    --reference $ref \
	    --outputs per_read_mods \
	    --output-directory ~/data/mdr/qc/megalodon/$samp/${samp} \
	    --write-mods-text \
	    --devices "cuda:0" \
	    --processes 36 > ~/data/mdr/qc/megalodon/$samp/${samp}/time.txt
fi

if [ $2 == megalodon_vanilla ] ; then
	mkdir -p ~/data/mdr/qc/megalodon/$samp/${samp}_vanilla

	time megalodon \
	    ~/data/mdr/qc/multiraw_sub/$samp \
	    --overwrite \
	    --guppy-server-path "/usr/bin/guppy_basecall_server" \
	    --reference $ref \
	    --outputs per_read_mods \
	    --output-directory ~/data/mdr/qc/megalodon/$samp/${samp}_vanilla \
	    --write-mods-text \
	    --devices "cuda:0" \
	    --processes 36 > ~/data/mdr/qc/megalodon/$samp/${samp}_vanilla/time.txt

fi


	       
    

