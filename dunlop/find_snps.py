import argparse
'''
take a vcf and find where the two non-ref seqs are different
take an annotation of the ref genome and find what genes are mutated
'''

def find_snps(vcffile):
    '''
    input vcf file path
    return list of lists [refchr, refposition, refbase, newbase, numreads_ref, numreads_alt]
    '''
    with open(vcffile) as f:
        content=f.read().split('\n')
        f.close()
    positions=[]
    for i in content:
        if len(i)>0:
            if i[0]!='#':
                info=i.split('\t')
                pile=info[9].split(':')
                ao=pile[5]
                ro=pile[3]
                positions.append([info[0], info[1], info[3], info[4], ro, ao])
    return(positions)


def prep_gff(gfffile):
    '''
    input gff
    return list: [chrname, start, end, genename, product, note]
    '''
    with open(gfffile) as f:
        content=f.read().split('\n')
        f.close()
    genes=[]
    for i in content:
        if len(i)>0:
            if i[0]!='#':
                info=i.split('\t')
                ##each 'gene' is redundant info since there's usually
                if info[3]!='region' and info[3]!='gene':
                    tags=info[8].split(';')
                    prodinfo=[x for x in tags if 'product=' in x]
                    geneinfo=[x for x in tags if 'gene=' in x]
                    noteinfo=[x for x in tags if 'Note=' in x]
                    if len(prodinfo)>0 and len(geneinfo)>0 and len(noteinfo)>0:
                        genes.append([info[0], info[3], info[4], ','.join(geneinfo), ','.join(prodinfo), ','.join(noteinfo)])
    return(genes)


def assign_snps(positions, genes):
    '''
    input snp positions and gene regions
    return snps in gene regions with gene info
    '''
    snpinfo=[]
    noncoding=0
    for snp in positions:
        oldlen=len(snpinfo)
        pos=int(snp[1])
        for j in genes:
            start=int(j[1])
            end=int(j[2])
            if snp[0]==j[0] and pos>start and pos<end:
                ##want [chr, pos, refbase, altbase, refreads, altreads, genename, geneprod, note]
                if len(snp[3].split(','))==1:
                    snpinfo.append([snp[0], snp[1], snp[2], snp[3], snp[4], snp[5], j[3], j[4], str(start), str(end)])
                else:
                    for alt in range(0,len(snp[3].split(','))):
                        altcount=snp[5].split(',')[alt]
                        altbase=snp[3].split(',')[alt]
                        snpinfo.append([snp[0], snp[1], snp[2], altbase, snp[4], altcount, j[3], j[4], str(start), str(end)])
        ##after iterating, check to see if a var has been added to the list. if not, then add it as an unknown
        if len(snpinfo)==oldlen:
            if len(snp[3].split(','))==1:
                snpinfo.append([snp[0], snp[1], snp[2], snp[3], snp[4], snp[5], 'unk', 'unk', 'unk', 'unk'])
            else:
                for alt in range(0,len(snp[3].split(','))):
                    altcount=snp[5].split(',')[alt]
                    altbase=snp[3].split(',')[alt]
                    snpinfo.append([snp[0], snp[1], snp[2], altbase, snp[4], altcount, 'unk', 'unk', 'unk', 'unk'])
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
    print(vcffile+','+str(total))

    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='annotate snps in a vcf using a gff')
    parser.add_argument('--vcf','-v',  help='vcf file made using like freebayes or something', type=str, required=True)
    parser.add_argument('--gff','-g',  help='annotation of the ref genome', type=str, required=True)
    parser.add_argument('--out','-o',  help='final csv file', type=str, required=True)
    args = parser.parse_args()

    main(args.gff, args.vcf, args.out)
