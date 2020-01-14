import argparse

def promoterSet(promoterfile, outfile):
    '''takes in the promotor set file and returns dict with promotor name and fa'''
    with open(promoterfile, 'r') as f:
        content=f.read().split('\n')
    f.close()
    fadict={}
    for i in content:
        if len(i)>0 and i[0]!='#':
            promotername=i.split('\t')[1].replace(' ', '_')
            fadict[promotername]=i.split('\t')[5]
    with open(outfile, 'w') as f:
        for i in fadict:
            f.write('>'+i+'\n')
            f.write(fadict[i]+'\n')
        f.close()
        
def pafToGff(paffile, outfile):
    '''takes paf file and promotor file and returns a gff'''
    with open(paffile, 'r') as f:
        content=f.read().split('\n')
    f.close()
    regions=[]
    ##[chr, database, feature, start, end, score, strand, phase[orf], attributes]
    regions.append(['chrom', 'database', 'feature', 'start', 'end', 'score', 'strand', 'orfphase', 'attributes'])
    for i in content[:-1]:
        name=i.split('\t')[0]
        chrom=i.split('\t')[5]
        start=int(i.split('\t')[7])
        end=int(i.split('\t')[8])
        score=i.split('\t')[11]
        strand=i.split('\t')[4]
        orfphase=i.split('\t')[9]
        attributes=';'.join(['gene=promoter','product='+name, 'Note=merged promoter from regulondb'])
        ##check if the region overlaps an existing region, merge the two into the biggest region
        changed=False
        for place, r in enumerate(regions):
            if chrom==r[0] and start<int(r[4]) and end>int(r[3]) and (start!=int(r[3]) or end!=int(r[4])):
                print(name+' overlaps')
                newstart=min(start, int(r[3]))
                newend=max(end, int(r[4]))
                prod=r[8].split(';')[1]
                newprod=prod.split('=')[1]+'_'+name
                newattributes=attributes=';'.join(['gene=promoter','product='+newprod , 'Note=merged promoter from regulondb'])
                regions[place]=[chrom, 'regulondb', 'promoter', str(newstart), str(newend), score, strand, orfphase, newattributes]
                changed=True
                break
        if changed==False: ##not the best to do it like this, but whatevs for now
            print(name+' does not overlap')
            attributes=';'.join(['gene=promoter','product='+name, 'Note=merged promoter from regulondb'])
            regions.append([chrom, 'regulondb', 'promoter', str(start), str(end), score, strand, orfphase, attributes])
    with open(outfile, 'w') as f:
        for r in regions:
            f.write('\t'.join(r)+'\n')
        f.close()

def main(infile, outfile):
    if outfile.endswith('.fa'):
        print('Converting promoter set to fasta')
        promoterSet(infile, outfile)
    elif outfile.endswith('.gff'):
        print('Converting paf file to gff')
        pafToGff(infile, outfile)
    else:
        print('outfile unclear')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='get promoter fasta from regulon promoter set, or get gff from paf of promotor alignments')
    parser.add_argument('--infile','-i',  help='promoter set file from regulon db or paf alignment', type=str, required=True)
    parser.add_argument('--outfile','-o',  help='has to end in .fa or .gff depending on what the expected output is', type=str, required=True)
    args = parser.parse_args()
    
    main(args.infile, args.outfile)
    
                                
