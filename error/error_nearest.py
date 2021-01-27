def tig_to_pos(mumsnps):
    '''get a dictionary of tigs that map to list of its positions of error'''
    with open(mumsnps, 'r') as f:
        content=f.read().splitlines()
    posdict={}
    for i in content:
        tig=i.split('\t')[10]
        if tig not in posdict:
            posdict[tig]=[int(i.split('\t')[0])]
        else:
            posdict[tig].append(int(i.split('\t')[0]))
    return posdict

def find_neighbor(posdict1, posdict2):
    '''Scan through list of positions and return list of distances
    now that I think about it, this maybe should be done by regions as well'''
    distances=[]
    for tig in posdict1:
        for i in posdict1[tig]:
            if tig in posdict2:
                diffs=[ abs(i-j) for j in posdict2[tig] ]
                distances.append(min(diffs))
            else:
                print(tig)
    return distances

def main(mumfile1, mumfile2, outfile):
    mum1=tig_to_pos(mumfile1)
    mum2=tig_to_pos(mumfile2)
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
    import argparse
    parser = argparse.ArgumentParser(description='get info about genomic distance to the nearest error')
    parser.add_argument('--mumfile1', '-m1', help='first mummer snp file', type=str)
    parser.add_argument('--mumfile2', '-m2', help='second mummer snp file', type=str)
    parser.add_argument('--outfile', '-o', help='output csv path', type=str)
    args = parser.parse_args()
    main(args.mumfile1, args.mumfile2,  args.outfile)
                        
