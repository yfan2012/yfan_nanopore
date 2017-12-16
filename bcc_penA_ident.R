datadir='~/Dropbox/Lab/carbapenem_r21/annotations/bucc/ident/'
refdir='/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs/'
writedir='~/Dropbox/Lab/carbapenem_r21/annotations/bucc/ident/'

tsvs=list.files(path=datadir,recursive=TRUE, pattern='*tsv')

positions=data.frame(gene=character(),
                     chromosome=character(),
                     identity=numeric(),
                     align_len=numeric(),
                     mismatch=numeric(),
                     gap_open=numeric(),
                     gene_start=numeric(),
                     gene_end=numeric(),
                     strain_start=numeric(),
                     strain_end=numeric(),
                     eval=numeric(),
                     bit=numeric(),
                     strain=character())

##for (i in c(15:25)) {
for (i in c(1:length(tsvs))) {
    strain=tsvs[i]
    strain_info=read.table(paste0(datadir,strain), stringsAsFactor=FALSE, sep='\t')
    strain_info$strain=strain
    colnames(strain_info)=c('gene', 'chr', 'identity', 'align_len', 'mismatch', 'gap_open', 'gene_start', 'gene_end', 'strain_start', 'strain_end', 'eval', 'bit', 'strain')
    ident=data.frame(gene=strain_info$gene, identity=strain_info$identity, alignment_length=strain_info$align_len)
    prefix=strsplit(strain, '[.]')[[1]][1]
    write.table(ident, file=paste0(writedir,prefix, '.csv'), sep=',', row.names=FALSE) 
}






