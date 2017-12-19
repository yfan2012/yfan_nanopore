from fasta_utils import fasta_dict, rev_comp

'''Take in blast positions and a reference genome, and write out a fasta'''

##I'm so tired. Just hardcoding info based on visual inspection of tsv. Srsly future yfan, pls forgive this sin.
##parC is reverse complement (blast results s.end is less than s.start)
gyrAloc=['>tig00000002_pilon', 3525263, 3527866]
parCloc=['>tig00000002_pilon', 1499655, 1501976]


refpath='/home/yfan/Dropbox/Lab/carbapenem_r21/annotations/fasta_corrected_assembly/illumina_corrected/170816_BUCC.pilon.fasta'

ref=fasta_dict(refpath)


gyrA=ref[gyrAloc[0]][gyrAloc[1]:gyrAloc[2]]
parCrev=ref[parCloc[0]][parCloc[1]:parCloc[2]]
parC=rev_comp(parCrev)


outdir='/home/yfan/Dropbox/Lab/carbapenem_r21/annotations/bucc/'


with open(outdir+'gyrA.fasta', 'w') as f:
    f.write('>tig2_gyrA\n')
    f.write(gyrA+'\n')

with open(outdir+'parC.fasta', 'w') as f:
    f.write('>tig2_parC\n')
    f.write(parC+'\n')

