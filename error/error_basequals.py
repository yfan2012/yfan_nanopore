from error_nearest import tig_to_pos
from error_venn_3samps import collapse
import pysam



def main(mumfile1, bamfile1, bamfile2, outfile):
    '''take in a mummer snps file and both bam files
    return quals from both bam files for each position listed in the mummer snp file
    it can be a downstream problem to compare the sets - making it easier this time too, by actually'''
    mum1=tig_to_pos(mumfile1)    
    sam1=pysam.AlignmentFile(bamfile1, 'rb')
    sam2=pysam.AlignmentFile(bamfile2, 'rb')
    qualinfo=[]
    for tig in mum1:
        print(tig)
        pos1=collapse(mum1[tig])
        for i in pos1:
            basepos=0
            quals=[]
            for q in sam1.pileup(tig, i[0], i[1], truncate=True):
                quals.append([tig, i[0], i[1], i[0]+basepos, sum(q.get_query_qualities())/len(q.get_query_qualities())])
                basepos+=1
            basepos=0
            for q in sam2.pileup(tig, i[0], i[1], truncate=True):
                quals[basepos].append(sum(q.get_query_qualities())/len(q.get_query_qualities()))
                basepos+=1
            qualinfo+=quals
    with open(outfile, 'w') as f:
        f.write(','.join(['tig', 'error_start', 'error_end', 'error_position', 'sampqual', 'otherqual'])+'\n')
        for i in qualinfo:
            f.write(','.join([i[0], str(i[1]), str(i[2]), str(i[3]), str(i[4]), str(i[5])])+'\n')
        f.close()

if __name__ == '__main__':
    import pysam
    import argparse
    parser = argparse.ArgumentParser(description='get basequalities of error blocks as defined by ')
    parser.add_argument('--mumfile', '-m', help='mummer snp file', type=str)
    parser.add_argument('--bamfile1', '-b1', help='bam file', type=str)
    parser.add_argument('--bamfile2', '-b2', help='second bamfile', type=str)
    parser.add_argument('--outfile', '-o', help='output csv path', type=str)
    args = parser.parse_args()
    main(args.mumfile, args.bamfile1, args.bamfile2, args.outfile)
                        
