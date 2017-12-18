import argparse
import re

parser=argparse.ArgumentParser(description='why do i code so slow')
parser.add_argument('--gff', '-g', type=str, required=True, help='full path to gff file')
args=parser.parse_args()

with open(args.gff) as f:
    content=f.readlines()


prefix=args.gff.split('/')[-1][:-4]
print prefix

extract=[]
for i in content:
    if i[0]=='#' :
        content.remove(i)
    else:
        rna_type=re.split('=|;', i)[1]
        if rna_type=='16S_rRNA':
            info=i.split('\t')
            extract.append([info[0], info[3], info[4]])

path=args.gff.split('/')[:-2]
path.append('bcc_refs/'+prefix+'.fa')
fasta='/'.join(path)

with open(fasta) as f:
    content=f.read().splitlines() 


dictfa={}
seq=''
name=''
for i in content:
    if i[0]=='>':
        dictfa[name]=seq
        name=i.split(' ')[0]
        seq=''
    else:
        seq+=i
dictfa[name]=seq


count=0
for i in extract:
    gene=dictfa['>'+i[0]][int(i[1]):int(i[2])]
    count_mod=str(count)
    genepath=args.gff.split('/')[:-2]
    genepath.append('bcc_16S_genes/'+prefix+count_mod+'.fa')
    genefa='/'.join(genepath)
    with open(genefa, 'w') as f:
        f.write('>'+i[0]+'\n')
        f.write(gene)
    count+=1
