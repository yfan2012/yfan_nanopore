def tig_to_pos(mumfile):
    '''get a dictionary of tigs that mapp to list of its positions of error'''
    with open(mumfile, 'r') as f:
        content=f.read().splitlines()
    posdict={}
    for i in content:
        tig=i.split('\t')[10]
        if tig not in posdict:
            posdict[tig]=[int(i.split('\t')[0])]
        else:
            posdict[tig].append(int(i.split('\t')[0]))
    return posdict

def get_bam_info(bamfile, posdict1,):
    '''take a bamfile and return basequal info'''
    samfile=pysam.AlignmentFile(bamfile, "rb")


def get_quals(tig, start, flank, baminfo):
    '''take a location on the reference and return list of all basequals near it'''
    
    


def main(mumfile, bamfile, outfile):
    mum=tig_to_pos(mumfile)


    print('snps in r9')
    dist1=find_neighbor(mum1, mum2)
    print('snps in r10')
    dist2=find_neighbor(mum2, mum1)
    with open(outfile, 'w') as f:
        f.write(','.join(['pore', 'distance'])+'\n')
        for i in dist1:
            f.write(','.join(['r9', str(i)])+'\n')
        for i in dist2:
            f.write(','.join(['r10',str(i)])+'\n')

if __name__ == '__main__':
    import pysam
    import argparse
    parser = argparse.ArgumentParser(description='get info about genomic distance to the nearest error')
    parser.add_argument('--mumfile', '-m', help='mummer snp file', type=str)
    parser.add_argument('--bamfile', '-b', help='bam file', type=str)
    parser.add_argument('--outfile', '-o', help='output csv path', type=str)
    args = parser.parse_args()
    main(args.mumfile, args.bamfile,  args.outfile)
                        
