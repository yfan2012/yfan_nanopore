import pysam

samfile=pysam.AlignmentFile("/atium/Data/Nanopore/cpowgs/bucc/170816_BUCC.sorted.bam", "rb")
tigs={'tig00000001':3156674, 'tig00000002':3822728, 'tig00000003':1198404, 'tig00000006':195378, 'tig00000010':49768, 'tig00000011':176624}

for i in tigs:
    total_nts=0
    for read in samfile.fetch(i, 1, tigs[i]):
        total_nts+=read.query_alignment_length
    cov=float(total_nts)/float(tigs[i])
    print i+'\t'+str(cov)
    
