import pysam
import sys
sys.path.insert(0, '/home/yfan/Code/utils')
from fasta_utils import fasta_dict
'''
mdr nt set of references contains most plasmids associated with the top40 species found by kraken
mdr ref set of references contains just reference level genomes
take just the plasmids of mdr nt, and the chrs from mdr ref
'''


mdr_nt='/atium/Data/ref/mdr_nt/mdr_nt.fa'
mdr_ref='/mithril/Data/Nanopore/projects/methbin/mdr/ref/mdr_refs.fa'

ntdict=fasta_dict(mdr_nt)
refdict=fasta_dict(mdr_ref)

ntred={}
rnameinfo={}

for i in ntdict:
    if 'plasmid' in i:
        ntred[i]=ntdict[i]
        acc=i.split(' ')[0][1:]
        rnameinfo[acc]='plasmid'

for i in refdict:
    if 'plasmid' not in i:
        ntred[i]=refdict[i]
        acc=i.split(' ')[0][1:]
        rnameinfo[acc]='chromosome'

        
outfile='/atium/Data/ref/mdr_nt/mdr_reduced.fa'
keyfile='/atium/Data/ref/mdr_nt/mdr_reduced_key.txt'

with open(outfile, 'w') as f:
    for i in ntred:
        f.write(i+'\n')
        f.write(ntred[i]+'\n')

with open(keyfile, 'w') as f:
    for i in rnameinfo:
        f.write('\t'.join([i, rnameinfo[i]])+'\n')
        
