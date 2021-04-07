import pysam
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='get transcriptome sequences that have been annotated by interproscan')
    parser.add_argument('--trinity','-t',  help='transcriptome', required=True, type=str)
    parser.add_argument('--interproscan','-i',  help='interproscan tsv', required=True, type=str)
    parser.add_argument('--outfile','-o',  help='output fasta', required=True, type=str)
    args = parser.parse_args()
    return args

def get_relevant_seqs(trindict, interlist):
    '''
    takes fastadict and interproscan list of seqnames
    returns fasta of only sequences included in the interproscan list
    '''
    needseqs={}
    for i in trindict:
        if i in interlist:
            needseqs[i]=trindict[i]
    return needseqs

def get_interproscan_names(intertsvfile):
    '''
    takes interproscan tsv file
    returns list of seqnames in that file
    '''
    namelist=[]
    with open(intertsvfile, 'r') as f:
        content=f.read().splitlines()
    for i in content:
        namelist.append(i.split('\t')[0].split('.')[0])
    interlist=set(namelist)
    return interlist

def main(trinfafile, intertsvfile, outfile):
    ##read in transcriptome fasta
    trin=pysam.FastaFile(trinfafile)
    xscripts=trin.references
    trindict={ x:trin.fetch(x) for x in xscripts }

    interlist=get_interproscan_names(intertsvfile)
    needseqs=get_relevant_seqs(trindict, interlist)

    with open(outfile, 'w') as f:
        for i in needseqs:
            f.write('>'+i+'\n')
            f.write(needseqs[i]+'\n')

if __name__ == '__main__':
    args=parseArgs()
    main(args.trinity, args.interproscan, args.outfile)

