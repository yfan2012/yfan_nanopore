##Read blast report of where genes should be
datadir='/atium/Data/Nanopore/cpowgs/genomes/var_profiles'

with open(datadir+'/ref_genes.tsv') as f:
    content=f.read().splitlines()

queries=[]
for i in content[1::6]:
    queries.append(i.split(' ')[2])

info=[]
for i in content[5::6]:
    info.append([i.split('\t')[8], i.split('\t')[9]])

ranges=dict(zip(queries, info))

ref='/mithril/Data/NGS/Reference/kpneumo/NC_016845.fasta'
with open(ref) as f:
    content=f.read().split('\n')


refgenes='/atium/Data/Nanopore/cpowgs/References/ref_gene_seqs'
for i in ranges:
    check=int(ranges[i][1])-int(ranges[i][0])
    if check > 0:
        seq=content[1][int(ranges[i][0]):int(ranges[i][1])]
    else:
        rev=content[1][int(ranges[i][1]):int(ranges[i][0])]
        seq=rev[::-1]
    with open(refgenes+'/'+i+'.fasta', 'w') as f:
        f.write('>'+i+' from HS11286\n')
        f.write(seq+'\n')
        
    
    
    
    
