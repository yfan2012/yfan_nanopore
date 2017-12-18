import os
import argparse

parser=argparse.ArgumentParser(description='get alignment percent')
parser.add_argument('--input', '-i', type=str, required=True, help='full path to flagstat.txt')
args=parser.parse_args()

with open(args.input) as f:
    content=f.readlines()

total=float(content[0].split(' ')[0])
aligned=float(content[4].split(' ')[0])

if aligned!=0:
    percent_align=aligned/total
else:
    percent_align=0
    
label=args.input.split('.')[0].split('/')[-1]

print ','.join([label, str(percent_align)])
