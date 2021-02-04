import argparse
import pysam
from Bio import SeqIO

def read_regions(regionfile):
    '''
    read in regions file as csv
    return list of regions [tig, start, end]
    '''
    with open(regionfile) as f:
        content=f.read().split('\n')
        f.close()
    regions=[]
    for i in content:
        if len(i)>0:
            info=i.split(',')
            regions.append([info[1], int(info[2]), int(info[3])])
    return regions

def filter_regions(regions, asmdict, end=0):
    '''
    filter list of lists regions by end 
    return list of lists that still should be breakpoints
    '''
    keepregions=[]
    for i in regions:
        ##end filtering: keep region for later deletion if any part is inside the end threshold
        lower=0+end
        upper=len(asmdict[i[0]])-end
        if i[2] > lower and i[1] <= upper:
            keepregions.append(i)
    return keepregions


def align_span(keepregions, npbam, minspan=5):
    '''
    take in list of lists regions and relevant bam
    return list of lists regions that are NOT spanned by any alignments
    '''
    for i in keepregions:
        lower=i[1]-5
        upper=i[2]+5
        

        

def break_tigs(asmfile, npbamfile, illbamfile, regionfile):
    '''
    take a fasta and list of zero coverage [tig, start, end]
    return seqio records of broken contigs
    option to disregard ends of contigs
    '''
    regions=read_regions(regionfile)
    asmdict=SeqIO.to_dict(SeqIO.parse(asmfile, 'fasta'))
    keepregions=filter_regions(regions, asmdict, end=100)

    npbam=pysam.AlignmentFile(npbamfile, 'rb')
    npregions=align_span(keepregions, npbam)
    
    illbam=pysam.AlignmentFile(illbamfile, 'rb')
    illregions=align_span(npregions, illbam)
    
    for i in tigs:
        
        
        


def main(asmfile, regionfile, outfile):

    newtigs=break_tigs(asmfile, regions)
    with open(outfile, 'w') as f:
        SeqIO.write(newtigs, f, 'fasta')
