#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/read_class/disco
datadir=/mithril/Data/Nanopore/projects/methbin/disco
ssddir=~/data/mdr/disco

#https://sra-download.ncbi.nlm.nih.gov/traces/sra71/SRZ/010032/SRR10032544/MinION_HP_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRZ/010032/SRR10032548/MinION_CP_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra31/SRZ/010032/SRR10032558/MinION_NG_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRZ/010032/SRR10032560/MinION_MH_NAT.tar.gz
#https://sra-download.ncbi.nlm.nih.gov/traces/sra25/SRZ/010032/SRR10032567/MinION_BA_NAT.tar.gz


if [ $1 == unpack ] ; then
    #samps='MinION_NG_NAT MinION_CP_NAT MinION_MH_NAT MinION_TP_NAT MinION_BF_NAT MinION_HP_NAT MinION_NG_NAT MinION_BA_NAT'
    samps=MinION_BA_NAT
    mkdir -p $ssddir/raw
    for i in $samps ;
    do
	tar -xzf $rawdir/$i.tar.gz -C $ssddir/raw
    done
fi




i=MinION_JM3O_NAT
if [ $1 == grab_metagenomes ] ; then
    aws s3 cp s3://sra-pub-src-1/SRR10139835/MinION_JM3T_NAT.tar.gz.1 $rawdir/
    aws s3 cp s3://sra-pub-src-2/SRR10038437/MinION_JM3O_NAT.tar.gz.1 $rawdir/
fi

if [ $1 == unpackmouse1 ] ; then
    mkdir -p $ssddir/raw
    tar -xzf $rawdir/$i.tar.gz.1 -C $ssddir/raw
fi

if [ $1 == ref_meta ] ; then
   mkdir -p $datadir/ref_meta
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/803/915/GCF_004803915.1_ASM480391v1/GCF_004803915.1_ASM480391v1_genomic.fna.gz -O $datadir/ref_meta/DUDU.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/688/845/GCF_001688845.2_ASM168884v2/GCF_001688845.2_ASM168884v2_genomic.fna.gz -O $datadir/ref_meta/MUIN.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/156/645/GCF_002156645.1_ASM215664v1/GCF_002156645.1_ASM215664v1_genomic.fna.gz -O $datadir/ref_meta/LAJO.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/803/695/GCF_004803695.1_ASM480369v1/GCF_004803695.1_ASM480369v1_genomic.fna.gz -O $datadir/ref_meta/MUGO.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/215/245/GCF_900215245.1_IMG-taxon_2617270901_annotated_assembly/GCF_900215245.1_IMG-taxon_2617270901_annotated_assembly_genomic.fna.gz -O $datadir/ref_meta/PSFL.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O $datadir/ref_meta/ESCO.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/412/675/GCF_000412675.1_ASM41267v1/GCF_000412675.1_ASM41267v1_genomic.fna.gz -O $datadir/ref_meta/PSPU.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz -O $datadir/ref_meta/SAEN.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/021/205/GCF_000021205.1_ASM2120v1/GCF_000021205.1_ASM2120v1_genomic.fna.gz -O $datadir/ref_meta/BACE.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/161/495/GCF_000161495.1_ASM16149v1/GCF_000161495.1_ASM16149v1_genomic.fna.gz -O $datadir/ref_meta/BATH.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/267/615/GCF_013267615.1_ASM1326761v1/GCF_013267615.1_ASM1326761v1_genomic.fna.gz -O $datadir/ref_meta/PADI.fa.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz -O $datadir/ref_meta/KLPN.fa.gz
fi

if [ $1 == make_ref_meta ] ; then
    cat $datadir/ref_meta/*.fa.gz > $datadir/ref_meta/meta10.fa.gz
    gunzip $datadir/ref_meta/meta10.fa.gz
fi

if [ $1 == move_mouse ] ; then
    mv $ssddir/megalodon/$i $datadir/megalodon/$i
    ##mv $ssddir/megalodon/${i}_meta $datadir/megalodon/${i}_meta
fi

