library(Biostrings)
library(BSgenome)

datadir='/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_vars/'
##refdir='/atium/Data/Nanopore/cpowgs/170816_BUCC/refs/'
refdir='/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_refs/'
writedir='/atium/Data/Nanopore/cpowgs/170816_BUCC/bcc_penA_genes/'

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
    positions=rbind(positions, strain_info)
}

write.table(positions, file='~/Dropbox/Lab/carbapenem_r21/annotations/bucc/bcc_ref_strains.csv', sep=',', row.names=FALSE)

positions=positions[positions$gene=='penA',]
positions$diff=positions$strain_end-positions$strain_start

strains=unique(positions$strain)


gene_seqs=character()
for (i in 1:dim(positions)[1]) {
    name=substr(positions$strain[i], 1, nchar(positions$strain[i])-15)    
    ref=readDNAStringSet(paste0(refdir, name, '.fa'), format='fasta') 
    first=min(positions$strain_start[i], positions$strain_end[i])
    last=max(positions$strain_start[i], positions$strain_end[i])

    newnames=as.character()
    for (j in 1:length(ref)) {
        newnames=c(newnames, strsplit(names(ref[j]), ' ')[[1]][1])
    }
    names(ref)=newnames
    
    wholechr=as.character(getSeq(ref, names=positions$chr[i]))
    seq=substr(wholechr, first, last)
    names(seq)=NULL

    if (positions$diff[i]<0) {
        chars=strsplit(seq, NULL)[[1]]
        seq=paste(rev(chars), collapse='')
    }

    gene_seqs=c(gene_seqs, seq)
}

positions$geneseqs=gene_seqs


for (i in strains) {
    aligns=sum(positions$strain==i)
    sname=substr(i, 1, nchar(i)-15)
    seqname=paste0('penA_', sname)

    ##This is a v janky way to check if the alignments overlap
    ranges=c()
    for (j in 1:aligns) {
        start=positions[positions$strain==i,]$gene_start[j]
        end=positions[positions$strain==i,]$gene_end[j]
        ranges=c(ranges, start:end)
    }
    cover=length(unique(ranges))
    overlap=length(ranges)-cover
    
    if (aligns==1) {
        newseq=DNAStringSet(positions$geneseqs[positions$strain==i])
        names(newseq)=seqname
        writeXStringSet(newseq, paste0(writedir, 'penA_', sname, '.fasta'),format="fasta", append=FALSE)
    } else if (overlap>500){
        subpositions=positions[positions$strain==i,]
        for (j in 1:dim(subpositions)[1]){
            newseq=DNAStringSet(subpositions$geneseqs[j])
            names(newseq)=seqname
            writeXStringSet(newseq, paste0(writedir, 'penA_', sname, as.character(j), '.fasta'),format="fasta", append=FALSE)
        }
    } else {
        totalseq=character()
        subpositions=positions[positions$strain==i,]
        order=sort(subpositions$gene_start)
        for (j in 1:length(order)) {
            totalseq=paste0(totalseq, subpositions$geneseqs[subpositions$gene_start==order[j]])
        }
        newseq=DNAStringSet(totalseq)
        names(newseq)=seqname
        writeXStringSet(newseq,paste0(writedir, 'penA_', sname, '.fasta'),format="fasta", append=FALSE)
    }
}    


