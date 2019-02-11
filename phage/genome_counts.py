'''get number of reads that aligned to each genome'''

def viral_counts(bampath):
    import pysam
    bamfile=pysam.AlignmentFile(bampath, 'r')
    xscript_counts={}
    for read in bamfile:
        ##if read.reference_id != -1 and read.reference_name in xscript_counts:
        if read.reference_name in xscript_counts:
            xscript_counts[read.reference_name] += 1
        else:
            xscript_counts[read.reference_name] = 1
    list_counts=[]
    for i in xscript_counts:
        list_counts.append([i, xscript_counts[i], float(xscript_counts[i])/float(sum(xscript_counts.values()))])
    list_counts.sort(key=lambda tup: tup[2], reverse=True)
    return list_counts


def main(bampath, refpath, outpath):
    counts=viral_counts(bampath)
    countsinfo=[[str(x[0]), str(x[1]), str(x[2])] for x in counts]
    
    with open(outpath, 'w') as f:
        f.write(','.join(['phage', 'counts', 'percent_of_aligned'])+ '\n')
        for i in countsinfo:
            f.write(','.join(i)+'\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='count how many reads align to each ref sequence')
    parser.add_argument('--bam','-i',  help='input bam path', type=str, required=True)
    parser.add_argument('--ref','-r',  help='input ref path (used to get seq names)', type=str, required=True)
    parser.add_argument('--csv','-o',  help='output csv path', type=str, required=True)
    args = parser.parse_args()
    
    main(args.bam, args.ref, args.csv)
