def fasta_dict(fastapath):
    '''fasta names are keys, sequences are values'''
    fasta={}
    with open(fastapath, 'r') as f:
        content=f.readlines()

    name=''
    seq=''
    for i in content:
        if i[0]=='>':
            if len(name)>0 and len(seq)>0:
                fasta[name]=seq
            name=i.strip('\n')
            seq=''
        else:
            seq+=i.strip('\n')
    fasta[name]=seq

    return(fasta)


def rev_comp(seq):
    key={'A': 'T', 'T':'A', 'G':'C',  'C':'G'}
    newseq=''
    for i in seq:
        newseq+=key[i]
    
    return newseq[::-1]
                            
