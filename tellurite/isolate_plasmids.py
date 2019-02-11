def list_tigs(blastfile):
    '''read the blast file and add all new contigs to the list'''
    with open(blastfile, 'r') as f:
        content=f.read().splitlines()
    hits=[]
    for i in content:
        if not i.startswith('#'):
            hits.append(i.split('\t')[1])
    tigs=list(set(hits))
    return tigs


def main(asmfile, blastfile, outfile):
    '''spit out ter plasmid only'''
    import sys
    sys.path.insert(0, '/home/yfan/Code/utils')
    from fasta_utils import fasta_dict
    contigs=fasta_dict(asmfile)
    tigs=['>'+x for x in list_tigs(blastfile)]
    with open(outfile, 'w') as f:
        for i in contigs:
            if i in tigs:
                f.write(i+'\n')
                f.write(contigs[i]+'\n')
                


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='get info about tellurite')
    parser.add_argument('--blast', '-b', help='blast outfmt 7 file', type=str)
    parser.add_argument('--asm', '-a', help='assembly', type=str)
    parser.add_argument('--outfile', '-o', help='output csv path', type=str)
    args = parser.parse_args()

    main(args.asm, args.blast, args.outfile)
