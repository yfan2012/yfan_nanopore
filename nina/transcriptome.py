import pysam
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description='get transcriptome sequences that have been annotated by interproscan')
    parser.add_argument('--trinity','-t',  help='transcriptome', required=True, type=str)
    parser.add_argument('--interproscan','-i',  help='interproscan tsv', required=True, type=str)
    parser.add_argument('--outfile','-o',  help='output fasta', required=True, type=str)
    parser.add_argument('--gofile','-g',  help='output tsv assigning GO terms to transcripts', required=True, type=str)
    parser.add_argument('--descriptfile','-d',  help='output tsv assigning pfam descriptions to transcripts', required=True, type=str)
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


def get_go_info(intertsvfile):
    '''
    takes interproscan tsv file
    GO info for each transcript
    use only ProSitePatterns
    '''
    with open(intertsvfile, 'r') as f:
        content=f.read().splitlines()
    goinfodict={}
    for i in content:
        info=i.split('\t')
        tname=info[0].split('_i')[0]
        if info[3]=='ProSitePatterns' and len(info)>13:
            if info[13]!='-':
                if tname in goinfodict:
                    goinfodict[tname].extend(info[13].split('|'))
                else:
                    goinfodict[tname]=info[13].split('|')
    for i in goinfodict:
        goinfodict[i]=list(set(goinfodict[i]))
    return goinfodict


def get_gen_info(intertsvfile):
    '''
    takes interproscan tsv file
    general descriptions from pfam
    use only ProSitePatterns
    '''
    with open(intertsvfile, 'r') as f:
        content=f.read().splitlines()
    geninfodict={}
    for i in content:
        info=i.split('\t')
        tname=info[0].split('_i')[0]
        if info[3]=='Pfam':
            if tname in geninfodict:
                geninfodict[tname].extend([info[5]])
            else:
                geninfodict[tname]=[info[5]]
    ##add mobidblite
    mobidict={}
    for i in content:
        info=i.split('\t')
        tname=info[0].split('_i')[0]
        if tname not in geninfodict:
            if tname in mobidict:
                mobidict[tname].append(info[5])
            else:
                mobidict[tname]=[info[5]]
    for i in mobidict:
        geninfodict[i]=mobidict[i]
    for i in geninfodict:
        geninfodict[i]=list(set(geninfodict[i]))
    return geninfodict


def main(trinfafile, intertsvfile, outfile, gofile, descriptfile):
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

    goinfodict=get_go_info(intertsvfile)
    with open(gofile, 'w') as f:
        for i in goinfodict:
            f.write('\t'.join([i, '|'.join(goinfodict[i])])+'\n')

    geninfodict=get_gen_info(intertsvfile)
    with open(descriptfile, 'w') as f:
        for i in geninfodict:
            f.write('\t'.join([i, '|'.join(geninfodict[i])])+'\n')

            
if __name__ == '__main__':
    args=parseArgs()
    main(args.trinity, args.interproscan, args.outfile, args.gofile, args.descriptfile)

