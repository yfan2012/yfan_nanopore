import glob
import re

coords=glob.glob("/home/yfan/Dropbox/Lab/carbapenem_r21/annotations/bucc/wga/*.tsv")


all_perc=[]
for tsv in coords:
    with open(tsv) as f:
        content=f.read().splitlines()

    org=re.split('/|\.', tsv)[-2]
        
    info=[]
    for i in content:
        info.append(i.split('\t'))

    del info[0:4]

    nts=0
    perc=[]
    for i in info:
        nts+=int(i[4])
        perc.append(float(i[6]))

    avgperc=sum(perc)/float(len(perc))

    cov=float(nts)/float(8300000)

    align_perc=cov*avgperc
    
    all_perc.append([org,str(align_perc), str(avgperc)])

with open('/home/yfan/Dropbox/Lab/carbapenem_r21/annotations/bucc/wga/est_align.csv', 'w') as f:
    for i in all_perc:
        f.write(','.join(x for x in i)+'\n')
