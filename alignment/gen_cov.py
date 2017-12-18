import argparse

parser=argparse.ArgumentParser(description='get genome coverage')
parser.add_argument('--cov', '-c', type=str, required=True, help='full path to cov file')
args=parser.parse_args()

with open(args.cov) as f:
    content=f.readlines()

info=[]
chrlengths=[]

for i in content:
    info.append(float(i.split('\t')[1])*float(i.split('\t')[2]))
    chrlengths.append(int(i.split('\t')[3]))

##Yes, doing chrlengths this way is awful, but I just need the answer rn

xcov=sum(info)/sum(set(chrlengths))
print set(chrlengths)
print '%f X coverage was acheived for this reference, according to alignments.' % xcov
