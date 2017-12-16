'''
take in fasta of penA gene, reverse complement, and write back to file
'''

import argparse
parser=argparse.ArgumentParser(description='fix my dumbness in R. Yes, I hate myself')
parser.add_argument('--pen', '-p', type=str, required=True, help='full path to penA fasta')
args=parser.parse_args()

pen='/atium/Data/Nanopore/cpowgs/170816_BUCC/penA_genes/penA_170816_BUCC.illumina.fasta'
with open(args.pen) as f:
    content=f.read().splitlines()

key={'A': 'T', 'T':'A', 'G':'C',  'C':'G'}

new=[]
for i in content:
    if i[0] != '>' :
        new.append(''.join([key[x] for x in i]))
new.insert(0, content[0])
        
with open(args.pen, 'w') as f:
    f.write('\n'.join(x for x in new))
