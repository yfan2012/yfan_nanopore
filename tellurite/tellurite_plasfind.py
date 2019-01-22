def plasmidlabel(plasmidfile):
    '''associate contigs with plasmids using abricate output of only plasmidfinder'''
    with open(plasmidfile, 'r') as f:
        content=f.read().splitlines()
    info=[x.split('\t') for x in content]
    chrplas={}
    ##in case one tig has multiple plasmid things on it
    done=[]
    for i in info[1::]:
        if i[1] in done:
            chrplas[i[1]]+=' and ' + i[4].replace('(', '').replace(')', '')
        else:
            chrplas[i[1]]=i[4].replace('(', '').replace(')', '')
        done.append(i[1])
    return chrplas

def readblast(blastfile):
    with open(blastfile, 'r') as f:
        content=f.read().splitlines()
    hits=[]
    for i in content:
        if not i.startswith('#'):
            hits.append(i.split('\t'))
    return hits

def main(blastfile, plasmidfile, genename, outfile, prefix):
    hits=readblast(blastfile)
    chrplas=plasmidlabel(plasmidfile)
    with open(outfile, 'w') as f:
        cols=['isolate','gene','accession','contig','plasmid','%ident','align_length','mismatch','gap_open', 'genestart','gene_end','contigstart','contigend','eval', 'bitscore']
        f.write(','.join(cols)+'\n')
        for i in hits:
            i.insert(0,genename)
            i.insert(0,prefix)
            if i[2] in chrplas:
                i.insert(3,chrplas[i[2]])
            else:
                i.insert(3,'unknown')
            f.write(','.join(i)+'\n')
    f.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='get info about tellurite')
    parser.add_argument('--blast', '-b', help='blast outfmt 7 file', type=str)
    parser.add_argument('--plasmid', '-p', help='plasmidfinder abricate file', type=str)
    parser.add_argument('--genename', '-g', help='genename label in blast file name', type=str)
    parser.add_argument('--prefix', '-n', help='prefix', type=str)
    parser.add_argument('--outfile', '-o', help='output csv path', type=str)
    args = parser.parse_args()
    main(args.blast, args.plasmid, args.genename, args.outfile, args.prefix)
