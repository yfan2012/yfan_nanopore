import argparse
'''
take gff
return how total number of coding bases
'''

def prep_gff(gfffile):
    '''
    input gff
    return list: [chrname, start, end, product]
    need exon info
    '''
    with open(gfffile) as f:
        content=f.read().split('\n')
        f.close()
    genes=[]
    for i in content:
        if len(i)>0:
            if i[0]!='#':
                info=i.split('\t')
                tags=info[8].split(';')
                ##ignore entries that are whole chromosomes, keep only the name and product
                if info[2] == 'CDS':
                    prodinfo=[x for x in tags if 'product=' in x]
                    parentinfo=[x for x in tags if 'Parent=' in x]
                    if len(prodinfo)>0 and len(parentinfo)>0:
                        genes.append([info[0], info[3], info[4], info[2]])
    return(genes)


def total_cds(genes):
    '''
    figure out how to add ranges
    take in prepped gff
    return number
    '''
    total={}
    for i in genes:
        if i[0] in total:
            total[i[0]]+=[*range(int(i[1]), int(i[2]))]
        else:
            total[i[0]]=[*range(int(i[1]), int(i[2]))]
        lengths=[]
    for i in total:
        lengths.append(len(set(total[i])))
    region_length=sum(lengths)
    return(region_length)


def main(gfffile):
    genes=prep_gff(gfffile)
    total=total_cds(genes)
    print(total)

    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='find reads that aligned exactly one time')
    parser.add_argument('--gff','-g',  help='annotation of the ref genome thrown into parsnp', type=str, required=True)
    args = parser.parse_args()

    main(args.gff)
