#!/bin/bash

ssddir=~/data/mdr/mdr
datadir=/mithril/Data/Nanopore/projects/methbin/mdr

if [ $1 == unpack ] ; then
    mkdir -p $ssddir/raw
    tardir=/uru/Data/Nanopore/projects/mdr

    tar -xzf $tardir/200708_mdr_stool16native.tar.gz -C $ssddir/raw
fi

    
if [ $1 == ref ] ; then
    mkdir -p $datadir/ref

    cd $datadir/ref
    ##bacteroides dorei (phocaeicola dorei)
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/387/545/GCF_902387545.1_UHGG_MGYG-HGUT-02478/GCF_902387545.1_UHGG_MGYG-HGUT-02478_genomic.fna.gz -O bacteroides_dorei.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/336/405/GCF_000336405.1_ASM33640v1/GCF_000336405.1_ASM33640v1_genomic.fna.gz -O enterococcus_faecium.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/131/755/GCF_014131755.1_ASM1413175v1/GCF_014131755.1_ASM1413175v1_genomic.fna.gz -O bacteroides_thetaiotamicron.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/407/245/GCF_000407245.1_Ente_aviu_ATCC14025_V2/GCF_000407245.1_Ente_aviu_ATCC14025_V2_genomic.fna.gz -O enterococcus_avium.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O escherichia_coli.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/286/525/GCF_001286525.1_BFBE1.1/GCF_001286525.1_BFBE1.1_genomic.fna.gz -O bacteroides_fragilis.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/825/GCF_000012825.1_ASM1282v1/GCF_000012825.1_ASM1282v1_genomic.fna.gz -O bacteroides_vulgatus.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz -O klebsiella_pneumoniae.fa.gz

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/949/455/GCF_000949455.1_ASM94945v1/GCF_000949455.1_ASM94945v1_genomic.fna.gz -O ruthenibacterium_lactatiformans.fa.gz
    
    cat *fa.gz > mdr_refs.fasta.gz
fi

    
