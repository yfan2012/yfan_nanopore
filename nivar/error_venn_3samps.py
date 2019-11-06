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

def collapse(poslist):
    '''Scan through list of positions and return list of ranges'''
    ##add pseudo location to the end so the while loop works easily
    ##if poslist[-1]-poslist[-2] == 1:
    poslist.append(999999999999999)
    diffs=[j-i for i,j in zip(poslist[:-1],poslist[1:])]
    ranges=[]
    i=0
    start=poslist[0]
    while i < len(diffs):
        if diffs[i] > 3:
            end=poslist[i]
            ranges.append([start-3, end+3])
            start=poslist[i+1]
            i+=1
        else:
            i+=1
    return ranges


def venn(posdict1, posdict2, posdict3):
    '''given two dictionaries of tigs and positions, return a list [tig, both, r1, r2]'''
    ##scan through positions and collapse anything that's less than two bases away
    vencounts=[]
    for tig in posdict1:
        if tig in posdict2 and tig in posdict3:
            collapsed1=collapse(posdict1[tig])
            collapsed2=collapse(posdict2[tig])
            collapsed3=collapse(posdict3[tig])
        allsamps1=0
        unique1=0
        only1_12=0
        only1_13=0
        for i in collapsed1:
            overlap2=[ i[0] <= j[1] and i[1] >= j[0] for j in collapsed2 ]
            overlap3=[ i[0] <= j[1] and i[1] >= j[0] for j in collapsed3 ]
            if sum(overlap2)>=1 and sum(overlap3)>=1:
                allsamps1+=1
            elif sum(overlap2)>=1 and sum(overlap3)<1:
                only1_12+=1
            elif sum(overlap3)>=1 and sum(overlap2)<1:
                only1_13+=1
            else:
                unique1+=1
        allsamps2=0
        unique2=0
        only2_21=0
        only2_23=0
        for i in collapsed2:
            overlap1=[ i[0] <= j[1] and i[1] >= j[0] for j in collapsed1 ]
            overlap3=[ i[0] <= j[1] and i[1] >= j[0] for j in collapsed3 ]
            if sum(overlap1)>=1 and sum(overlap3)>=1:
                allsamps2+=1
            elif sum(overlap1)>=1 and sum(overlap3)<1:
                only2_21+=1
            elif sum(overlap3)>=1 and sum(overlap1)<1:
                only2_23+=1
            else:
                unique2+=1
        allsamps3=0
        unique3=0
        only3_31=0
        only3_32=0
        for i in collapsed3:
            overlap1=[ i[0] <= j[1] and i[1] >= j[0] for j in collapsed1 ]
            overlap2=[ i[0] <= j[1] and i[1] >= j[0] for j in collapsed2 ]
            if sum(overlap1)>=1 and sum(overlap2)>=1:
                allsamps3+=1
            elif sum(overlap1)>=1 and sum(overlap2)<1:
                only3_31+=1
            elif sum(overlap2)>=1 and sum(overlap1)<1:
                only3_32+=1
            else:
                unique3+=1
        ##if two distinct regions in one map to the same region in the other, the intersection of the venn diagram can e different
        ##force them to be the smaller - this essentially combines the two distinct regions
        allsamps=min([allsamps1, allsamps2, allsamps3])
        only12=min([only1_12, only2_21])
        only13=min([only1_13, only3_31])
        only23=min([only2_23, only3_32])
        vencounts.append([tig, str(allsamps), str(only12), str(only13), str(only23), str(unique1), str(unique2), str(unique3)])
    return vencounts

def main(mumfile1, mumfile2, mumfile3, outfile):
    mum1=tig_to_pos(mumfile1)
    mum2=tig_to_pos(mumfile2)
    mum3=tig_to_pos(mumfile3)
    counts=venn(mum1, mum2, mum3)
    with open(outfile, 'w') as f:
        f.write(','.join(['tig', 'all', 'only12', 'only13', 'only23','unique1','unique2','unique3'])+'\n')
        for i in counts:
            f.write(','.join(i)+'\n')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='get info about which mummer snps are common between two queries aligned to the same ref')
    parser.add_argument('--mumfile1', '-m1', help='first mummer snp file', type=str)
    parser.add_argument('--mumfile2', '-m2', help='second mummer snp file', type=str)
    parser.add_argument('--mumfile3', '-m3', help='third mummer snp file', type=str)
    parser.add_argument('--outfile', '-o', help='output csv path', type=str)
    args = parser.parse_args()
    main(args.mumfile1, args.mumfile2, args.mumfile3,  args.outfile)
                        
