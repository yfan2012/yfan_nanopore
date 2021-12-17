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
    
    cat *fa.gz > mdr_refs.fa.gz
fi

dbdir=/atium/Data/ref/bacteria
if [ $1 == get_refseq_info ] ; then
    wget -O $dbdir/assembly_summary.txt ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
fi

if [ $1 == filt_refseq_info ] ; then
    awk -F '\t' '{if($12=="Complete Genome") print $20}' $dbdir/assembly_summary.txt > $dbdir/assembly_summary_complete_genomes.txt
    awk -F '\t' '{if($5!="na") print $20}' $dbdir/assembly_summary.txt > $dbdir/assembly_summary_rep_genomes.txt
fi


if [ $1 == dl_bacteria ] ; then
    mkdir -p $dbdir/genomes
    for next in $(cat $dbdir/assembly_summary_rep_genomes.txt) ;
    do
	name=`echo $next | rev | cut -d / -f 1 | rev `
	wget -P $dbdir/genomes "$next"/${name}_genomic.fna.gz ;
    done
fi

if [ $1 == cat_bacteria ] ; then
    ##some genomes corrupted? had to go thru and figure out which
    cat $dbdir/genomes/*fna.gz > $dbdir/all_bacteria_refs.fa.gz
    gunzip $dbdir/all_bacteria_refs.fa.gz
fi

if [ $1 == makeblastdb ] ; then
    mkdir -p $dbdir/blastdb
    makeblastdb \
	-in $dbdir/all_bacteria_refs.fa \
	-out $dbdir/blastdb/all_bacteria_refs \
	-dbtype nucl
fi

if [ $1 == get_fa_headers ] ; then
    ##species name is encoded in header most of the time
    grep '>' $dbdir/all_bacteria_refs.fa > $dbdir/all_bacteria_refs_faheaders.txt
    sed -i -e 's/>//g' $dbdir/all_bacteria_refs_faheaders.txt
    sed -i -e 's/[][]//g' $dbdir/all_bacteria_refs_faheaders.txt
fi


if [ $1 == get_tax_info ] ; then
    mkdir -p $dbdir/taxonomy
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
	 -O $dbdir/taxonomy/taxdump.tar.gz
    tar -xzf $dbdir/taxonomy/taxdump.tar.gz \
	-C $dbdir/taxonomy
fi

if [ $1 == clean_taxo ] ; then
    sed -e 's/[\t]*//g' $dbdir/taxonomy/names.dmp |\
	sed -e 's/|/,/g' | \
	sed -e 's/"//g' > $dbdir/taxonomy/names_clean.dmp

    sed -e 's/[\t*]//g' $dbdir/taxonomy/nodes.dmp |\
	sed -e 's/|/,/g' > $dbdir/taxonomy/nodes_clean.dmp
fi

prefix=200708_mdr_stool16native
if [ $1 == getplasmids ] ; then
    acclist=`awk '(NR>1) {print $13}' $datadir/amr/$prefix.plasmidfinder.tsv`
    for i in $acclist ;
    do
	esearch -db nuccore -query $i | efetch -format fasta > $datadir/ref/$i.fa
    done
fi

if [ $1 == lastfaecalis ] ; then
    ##one of the plasmid finder hits above I think is this facalis
    for i in 0 1 2 3 4 5 6 ;
    do
	acc=CP00662${i}
	esearch -db nuccore -query $acc | efetch -format fasta > $datadir/ref/$acc.fa
    done
fi

if [ $1 == mergerefs ] ; then
    cat $datadir/ref/*.fa > $datadir/ref/mdr_refs_plas.fa
fi

if [ $1 == mummer ] ; then
    mkdir -p $datadir/ref/mummer_self

    nucmer \
	-p $datadir/ref/mummer_self/mdr_refs_plas \
	$datadir/ref/mdr_refs_plas.fa \
	$datadir/ref/mdr_refs_plas.fa

    mummerplot --postscript \
	       -p $datadir/ref/mummer_self/mdr_refs_plas \
	       $datadir/ref/mummer_self/mdr_refs_plas.delta

    dnadiff \
	-p $datadir/ref/mummer_self/mdr_refs_plas \
	$datadir/ref/mdr_refs_plas.fa \
	$datadir/ref/mdr_refs_plas.fa

    show-coords -l -T $datadir/ref/mummer_self/mdr_refs_plas.delta \
		> $datadir/ref/mummer_self/mdr_refs_plas.sc.tsv
fi


if [ $1 == dl_nt_ncbi ] ; then
    ##download nt database for blasting
    ##written after the fact - downloaded on board winry first and moved, but this is the line of code
    mkdir -p /atium/Data/ref/ncbi
    cd /atium/Data/ref/ncbi

    update_blastdb.pl --decompress --num_threads 12 nt
fi


if [ $1 == dl_nt_fa ] ; then
    ##download the fasta of the nt database
    ##done at command line and recorded after the fact
    wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz -O /atium/Data/ref/ncbi/nt.gz
fi


if [ $1 == minimap_idx ] ; then
    minimap2 -t 36 -d /atium/Data/ref/ncbi/nt.mmi /atium/Data/ref/ncbi/nt.gz
fi

dbdir=/atium/Data/ref

if [ $1 == get_mdr_nt ] ; then
    ##try a different method to get all bacteria
    ##ecoli taxid is 562
    ##kleb taxid is 573
    mkdir -p $dbdir/mdr_nt

    ##get taxid from kraken
    awk '$5 != 9606 {print $5}' $datadir/kraken/$prefix.report.top40.txt \
	> $dbdir/mdr_nt/mdr_nt.taxid.txt

    ##get accessions
    touch $dbdir/mdr_nt/mdr_nt.acc.txt
    part=`echo $i | rev | cut -d _ -f 1 | rev`
    while read p; do
	echo $p
	blastdbcmd -db $dbdir/ncbi/nt -outfmt "%T %t %a" -entry all \
	    | awk -v var="$p" '$1 == var {print $0}' \
	    | grep 'plasmid\|chromosome\|genome' \
	    | grep complete \
	    | grep -v -w gene \
	    | grep -v -w genes \
	    | awk '{print $NF}' >> $dbdir/mdr_nt/mdr_nt.acc.txt
    done < $dbdir/mdr_nt/mdr_nt.taxid.txt
fi


if [ $1 == dl_fasta_nt ] ; then
    ##actually download sequences
    blastdbcmd \
	-entry_batch $dbdir/mdr_nt/mdr_nt.acc.txt \
	-db $dbdir/ncbi/nt \
	-out $dbdir/mdr_nt/mdr_nt.fa
fi


if [ $1 == CATdb ] ; then
    dbdir=$dbdir/CATdb
    mkdir -p $dbdir
    mkdir -p $dbdir/db
    mkdir -p $dbdir/tax
    
    CAT prepare --fresh \
	-d $dbdir/db \
	-t $dbdir/tax
    
    ##running CAT in a preprocess.sh
fi


    
if [ $1 == hiC_copy ] ; then
    mkdir -p $datadir/hiC/
    mkdir -p $datadir/hiC/clusters

    cp ~/gdrive/mdr/mdr/181127_stool_illumina_binning_rerun/proxiwrap_clusters/*fasta $datadir/hiC/clusters/
fi
