#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/disco
ssddir=~/data/mdr/disco

#https://sra-download.ncbi.nlm.nih.gov/traces/sra71/SRZ/010032/SRR10032544/MinION_HP_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRZ/010032/SRR10032548/MinION_CP_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra31/SRZ/010032/SRR10032558/MinION_NG_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRZ/010032/SRR10032560/MinION_MH_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRZ/010032/SRR10032567/MinION_BA_NAT.tar.gz

#samps='MinION_NG_NAT MinION_CP_NAT MinION_MH_NAT MinION_TP_NAT MinION_BF_NAT MinION_HP_NAT MinION_NG_NAT MinION_BA_NAT'
samps=MinION_BA_NAT

if [ $1 == unpack ] ; then
    for i in $samps ;
    do
	tar -xzf $rawdir/$i.tar.gz -C $ssddir
    done
fi


