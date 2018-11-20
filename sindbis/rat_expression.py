'''Get rat expression from transcriptome alignments'''

def xscript_counts(bampath):
    import pysam
    bamfile=pysam.AlignmentFile(bampath, 'r')
    xscript_counts={}
    for read in bamfile:
        if read.reference_name in xscript_counts:
            xscript_counts[read.reference_name] += 1
        else:
            xscript_counts[read.reference_name] = 1
    list_counts=[]
    for i in xscript_counts:
        list_counts.append([i, xscript_counts[i], float(xscript_counts[i])/float(sum(xscript_counts.values()))])
    list_counts.sort(key=lambda tup: tup[2], reverse=True)
    return list_counts

def seq_functions(refpath):
    seqfunc={}
    seqfunc[None]='unaligned'
    with open(refpath, 'r') as f:
        for line in f:
            if line[0] == '>':
                words=line.split(' ')
                seqfunc[words[0][1:]]=' '.join(words[1:]).rstrip()
    return seqfunc


def main(bampath, refpath, outpath):
    seqfunc=seq_functions(refpath)
    counts=xscript_counts(bampath)

    ##there should be a better way to do this:
    countsinfo=[[str(x[0]), seqfunc[x[0]], str(x[1]), str(x[2])] for x in counts]
        
    
    with open(outpath, 'w') as f:
        f.write(','.join(['seq_accession', 'function', 'counts', 'percent_of_aligned'])+ '\n')
        for i in countsinfo:
            f.write(','.join(i)+'\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='count how many drna reads align to each transcript in the rat xscriptome')
    parser.add_argument('--bam','-i',  help='input bam path', type=str, required=True)
    parser.add_argument('--ref','-r',  help='input ref path (used to get transcript functions)', type=str, required=True)
    parser.add_argument('--csv','-o',  help='output csv path', type=str, required=True)
    args = parser.parse_args()
    
    main(args.bam, args.ref, args.csv)
