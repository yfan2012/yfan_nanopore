#!/bin/bash

datadir=~/Dropbox/yfan/nina_fungus/assemblies
tmp=~/tmp/assemblies

if [ $1 == mum_assemblers ] ; then
    ##mummer assemblies against each other. Try to ditch wtdbg2 junk contigs in some way
    mkdir -p ~/mummer_assemblers

    ##clean up contig names in fastas and copy to a place where mummer can find it
    sed -i -e 's/_pilon//g' $datadir/*pilon/*fasta
    rm -r ~/tmp
    mkdir ~/tmp
    cp -r $datadir ~/tmp/
    
    for i in st31 st90853 ;
    do
	##mummer is being stupid about file paths. so write results to home first, then move the folder when everything is done
	nucmer -p ~/mummer_assemblers/$i $tmp/wtdbg2_pilon/$i*.fasta $tmp/canu_pilon/$i*20.fasta

	mummerplot --filter --fat --png -p ~/mummer_assemblers/$i.layout ~/mummer_assemblers/$i.delta -R $tmp/wtdbg2_pilon/$i*.fasta -Q $tmp/canu_pilon/$i*20.fasta
	mummerplot --filter --fat --png -p ~/mummer_assemblers/$i ~/mummer_assemblers/$i.delta

	dnadiff -p ~/mummer_assemblers/$i $tmp/wtdbg2_pilon/$i*.fasta $tmp/canu_pilon/$i*20.fasta
    done
    rm -rf $datadir/mummer_assemblers
    mv ~/mummer_assemblers $datadir/
fi

    
