'''
Quick check of penA gene coverage on pilon corrected assembly
'''

import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

samfile=pysam.AlignmentFile("/atium/Data/Nanopore/cpowgs/170816_BUCC/pilon/btbam/170816_BUCC.sorted.bam", "rb")

penAcov=[]
penApos=[]
for pileupcolumn in samfile.pileup("tig00000001", 2400580, 2402118):
    penAcov.append(pileupcolumn.n)
    penApos.append(pileupcolumn.pos)
    
f=plt.figure()
plt.scatter(penApos[200:len(penApos)-200], penAcov[200:len(penAcov)-200], s=.5)
plt.ylabel('Coverage', fontsize=8)
plt.xlabel('Contig/penA Position', fontsize=8)
plt.ylim(0, 200)
plt.title('penA coverage of illumina reads on hybrid assembly')
plt.tick_params(labelsize=5)

pp=PdfPages('/home/yfan/Dropbox/Lab/carbapenem_r21/annotations/etc/170816_BUCC_penA_pilon_cov.pdf')
pp.savefig(f)
pp.close()


