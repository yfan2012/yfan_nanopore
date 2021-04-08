import pysam

'''
quick check for unique genes at each step
'''


datadir='/pym/Data/Nanopore/projects/prolificans/rna'
intertsvfile=datadir+'/interproscan/Trinity.fasta.transdecoder_clean.pep.tsv'
tdecodefile=datadir+'/transdecoder/Trinity.fasta.transdecoder.pep'
trinfile=datadir+'/trinity/Trinity.fasta'

##interproscan
intergenes=[]
with open(intertsvfile, 'r') as f:
    content=f.read().splitlines()
    
for i in content:
    intergenes.append(i.split('_i')[0])
    
interunique=list(set(intergenes))


##tdecode
tdecodefa=pysam.FastaFile(tdecodefile)
tdecodegenes=[]
for i in tdecodefa.references:
    tdecodegenes.append(i.split('_i')[0])
    
tdecodeunique=list(set(tdecodegenes))


##trinity
trinfa=pysam.FastaFile(trinfile)
tringenes=[]
for i in trinfa.references:
    tringenes.append(i.split('_i')[0])

trinunique=list(set(tringenes))

print('Trinity:'+str(len(trinunique)))
print('transdecoder:'+str(len(tdecodeunique)))
print('interproscan:'+str(len(interunique)))


                              
