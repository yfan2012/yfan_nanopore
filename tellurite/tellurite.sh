#!/bin/bash

datadir=/kyber/Data/Nanopore/projects/tellurite
dbox=~/Dropbox/Timplab_Data/tellurite
abriout=$dbox/abricate
blastout=$dbox/blast_ter
repout=$dbox/report
prokout=$dbox/prokka

if [ $1 == get_ter_genes ] ; then
    ##get all ter c and d seqs using trick from geo
    esearch -db protein -query 'txid573[orgn] AND terC[prot]' | efetch -format fasta_cds_na > $datadir/ref/cds.terC.fa
    esearch -db protein -query 'txid573[orgn] AND terD[prot]' | efetch -format fasta_cds_na > $datadir/ref/cds.terD.fa
fi


if [ $1 == plasfind ] ; then
    for qual in raw pilon ;
    do
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    abricate --db plasmidfinder --quiet $i > $abriout/$qual/$prefix.plasmidfinder.tsv
	done
    done
fi

if [ $1 == virfind ] ; then
    for i in $datadir/assemblies/pilon/*fasta ;
    do
	prefix=`basename $i .fasta`
	abricate --db vfdb --quiet $i > $abriout/pilon/$prefix.vfdb.tsv
    done
fi

if [ $1 == mkblastdb ] ; then
    for qual in raw pilon ;
    do
	mkdir -p $datadir/blast_db/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    makeblastdb -in $i -out $datadir/blast_db/$qual/$prefix.db -dbtype nucl
	done
    done
fi

if [ $1 == search_ter ] ; then
    for qual in raw pilon ;
    do
	mkdir -p $blastout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
    	    prefix=`basename $i .fasta`
	    blastn -query $datadir/ref/cds.terC.fa -db $datadir/blast_db/$qual/$prefix.db -outfmt 7 -out $blastout/$qual/$prefix.terC.tsv
	    sed -i -e 's/_pilon//g' $blastout/$qual/$prefix.terC.tsv
	    blastn -query $datadir/ref/cds.terD.fa -db $datadir/blast_db/$qual/$prefix.db -outfmt 7 -out $blastout/$qual/$prefix.terD.tsv
	    sed -i -e 's/_pilon//g' $blastout/$qual/$prefix.terD.tsv
	done
    done
fi

if [ $1 == report ] ; then
    for qual in raw pilon ;
    do
	mkdir -p $repout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    echo $prefix
	    sed -i -e 's/_pilon//g' $abriout/$qual/$prefix.plasmidfinder.tsv
	    python ~/Code/yfan_nanopore/tellurite/tellurite_plasfind.py -b $blastout/$qual/$prefix.terD.tsv -p $abriout/$qual/$prefix.plasmidfinder.tsv -g terD -n $prefix -o $repout/$qual/$prefix.terD.report.csv
	    python ~/Code/yfan_nanopore/tellurite/tellurite_plasfind.py -b $blastout/$qual/$prefix.terC.tsv -p $abriout/$qual/$prefix.plasmidfinder.tsv -g terC -n $prefix -o $repout/$qual/$prefix.terC.report.csv
	done
    done
fi

if [ $1 == report_vfdb ] ; then
    for qual in pilon ;
    do
	mkdir -p $repout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    echo $prefix
	    sed -i -e 's/_pilon//g' $abriout/$qual/$prefix.plasmidfinder.tsv
	    python ~/Code/carbapenem_r21/report/abricate_report.py -c $abriout/$qual/$prefix.vfdb.tsv -p $abriout/$qual/$prefix.plasmidfinder.tsv -o $repout/$qual/$prefix.vfdbfind.report.csv
	done
    done
fi


if [ $1 == prokka ] ; then
    ##for qual in raw pilon ;
    for qual in pilon ;
    do
	mkdir -p $prokout/$qual
	for i in $datadir/assemblies/$qual/*fasta ;
	do
	    prefix=`basename $i .fasta`
	    prokka --outdir $prokout/$qual --genus Klebsiella --usegenus --prefix $prefix --cpus 12 --force $i &> $prokout/$qual/$prefix.log.txt
	done
    done
fi



if [ $1 == isolate_plas ] ; then
    for gene in terC terD ;
    do
	for i in /home/yfan/Dropbox/Timplab_Data/tellurite/assemblies/pilon/*fasta ; 
	do
	    prefix=`basename $i .fasta`
	    python ~/Code/yfan_nanopore/tellurite/isolate_plasmids.py -b $blastout/pilon/$prefix.$gene.tsv -a $i -o $dbox/assemblies/ter_plasmid_pilon/$prefix.$gene.tigs.fasta
	done
    done
fi


parsdir=$dbox/parsnp_plas
if [ $1 == parsnp ] ; then
    parsnp -r $dbox/assemblies/ter_plasmid_pilon/KLPN_2156.pilon.10.tigs.fasta -d $dbox/assemblies/ter_plasmid_pilon -p 12 -o $parsdir/parsnp1
    harvesttools -i $parsdir/parsnp1/parsnp.ggr -V $parsdir/parsnp1/parsnp.vcf

    parsnp -r $dbox/assemblies/ter_plasmid_pilon/KLPN_1165.pilon.10.tigs.fasta -d $dbox/assemblies/ter_plasmid_pilon -p 12 -o $parsdir/parsnp2
    harvesttools -i $parsdir/parsnp2/parsnp.ggr -V $parsdir/parsnp2/parsnp.vcf

    parsnp -r $dbox/assemblies/ter_plasmid_pilon/KLPN_1469.pilon.10.tigs.fasta -d $dbox/assemblies/ter_plasmid_pilon -p 12 -o $parsdir/parsnp3
    harvesttools -i $parsdir/parsnp3/parsnp.ggr -V $parsdir/parsnp3/parsnp.vcf

fi


if [ $1 == prokka_shorten ] ; then
    for i in $prokout/pilon/*.gff ;
    do
	prefix=`basename $i .gff`
	tig=`grep terC $i | cut -d $'\t' -f 1`
	if [ ! -z "$tig" ] ; then
	   grep 'gene=' $i | grep $tig > $prokout/pilon/$prefix.short.gff
	fi
	echo $prefix
	echo $tig
    done
fi

if [ $1 == prokka_dumb ] ; then
    for i in $prokout/pilon/*.short.gff ;
    do
	prefix=`basename $i .short.gff`
	if [ -f $prokout/pilon/$prefix.genelist.txt ] ; then
	    rm $prokout/pilon/$prefix.genelist.txt
	fi
	touch $prokout/pilon/$prefix.genelist.txt
	while read tig prod cds start end dot strand zero info ;
	do
	    gene=`echo $info | cut -d ';' -f 3 | cut -d '=' -f 2`
	    echo $gene >> $prokout/pilon/$prefix.genelist.txt
	done < "$i"
    done
fi

