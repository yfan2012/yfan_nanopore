import argparse
'''
take a parsnp vcf and find where the two non-ref seqs are different
take an annotation of the ref genome and find what genes are mutated
'''

def find_snps(vcffile):
    '''
    input vcf file path
    return list of lists [refchr, refposition, refbase, birdbase, ptbase]
    only includes positions where the two samp genomes aren't the same
    '''
    with open(vcffile) as f:
        content=f.read().split('\n')
        f.close()
    positions=[]
    for i in content:
        if len(i)>0 :
            if i[0]!='#':
                info=i.split('\t')
                if info[10]!=info[11]:
                    bases=info[4].split(',')
                    bases.insert(0, info[3])
                    positions.append([info[0], info[1], info[3], bases[int(info[10])],bases[int(info[11])]])
    return(positions)


def prep_gff(gfffile):
    '''
    input gff
    return list: [chrname, start, end, product]
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
                if not any('chromosome' in x for x in tags):
                    prodinfo=[x for x in tags if 'product=' in x]
                    parentinfo=[x for x in tags if 'Parent=' in x]
                    if len(prodinfo)>0 and len(parentinfo)>0:
                        genes.append([info[0], info[3], info[4], ','.join(parentinfo), ','.join(prodinfo)])
    return(genes)


def assign_snps(positions, genes):
    '''
    input snp positions and gene regions
    return snps in gene regions, with gene info
    '''
    snpinfo=[]
    noncoding=0
    for snp in positions:
        ##h8 this nested for lop nonsense, but i am too cranky to think of a better way rn
        pos=int(snp[1])
        for j in genes:
            start=int(j[1])
            end=int(j[2])
            if snp[0]==j[0] and pos>start and pos<end:
                ##want [chr, pos, refbase, birdbase, ptbase, geneparent, geneprod]
                snpinfo.append([snp[0], snp[1], snp[2], snp[3], snp[4], j[3], j[4]])
    return(snpinfo)


def main(gfffile, vcffile, outfile):
    positions=find_snps(vcffile)
    genes=prep_gff(gfffile)
    total=len(positions)
    a=assign_snps(positions, genes)
    with open(outfile, 'w') as f:
        for i in a:
            f.write(','.join(i)+'\n')
        f.close()
    print('Total of '+str(total)+' snps')
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='find reads that aligned exactly one time')
    parser.add_argument('--vcf','-v',  help='parsnp vcf file', type=str, required=True)
    parser.add_argument('--gff','-g',  help='annotation of the ref genome thrown into parsnp', type=str, required=True)
    parser.add_argument('--out','-o',  help='final csv file', type=str, required=True)
    args = parser.parse_args()

    main(args.gff, args.vcf, args.out)
