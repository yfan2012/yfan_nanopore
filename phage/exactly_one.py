def exactly_one_align(bampath, newbam):
    '''Find reads that have exactly one alignment by minimap2'''
    import pysam
    bamfile=pysam.AlignmentFile(bampath, 'r')
    not_unique=[]

    numalign=0
    for i in bamfile:
        if i.is_secondary == True:
            not_unique.append(i.query_name)
        numalign+=1
        
    bamfile=pysam.AlignmentFile(bampath, 'r')        
    uniquealignments=pysam.AlignmentFile(newbam, 'wb', template=bamfile)
    numunique=0
    for i in bamfile:
        if i.query_name not in not_unique:
            uniquealignments.write(i)
            numunique+=1
    uniquealignments.close()

    return [numalign, numunique]



def main(bampath, newbam):
    '''get reads that aligned exactly 1 time'''
    aligned=exactly_one_align(bampath, newbam)
    print(aligned[0])
    print(aligned[1])
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='find reads that aligned exactly one time')
    parser.add_argument('--bamfile','-i',  help='input bam path', type=str, required=True)
    parser.add_argument('--newfile','-o',  help='output bam path', type=str, required=True)
    args = parser.parse_args()
    
    main(args.bamfile, args.newfile)
