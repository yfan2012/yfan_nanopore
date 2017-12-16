library(argparse)

parser=ArgumentParser()
parser$add_argument('-d', '--directory', type='character', required=TRUE, help='isolate dir')
parser$add_argument('-i', '--illumina', type='character', required=TRUE, help='full prefix to illumina fq files ending in isolate number')
parser$add_argument('-c', '--canupre', type='character', required=TRUE, help='canu prefix')
args=parser$parse_args()

pilondir=file.path(args$directory,'pilon')
dir.create(file.path(pilondir))
outdir=args$output

##make index for alignment to canu contigs

dir.create(file.path(pilondir, "btidx"))
dir.create(file.path(pilondir, "btbam"))

##intially corrected is equal to canu assembly
system(paste0('cp ', file.path(args$directory, 'canu_assembly', paste0(args$canupre, '.contigs.fasta')), ' ', pilondir, '/btidx/'))
genome=file.path(pilondir, 'btidx', paste0(args$canupre, '.contigs.fasta'))


for (i in 1:10) {

    print('=================================================================================================================================================================')
    print(paste0('On round ', as.character(i)))

    ##build bowtie index and send to index dir made above, align and send to dir made above
    system(paste0("bowtie2-build ", genome, " ", pilondir, '/btidx/', args$canupre))
    system(paste0("bowtie2 -p 12 -x ", pilondir, "/btidx/", args$canupre, " -1 ", args$illumina, "_R1_001.fastq.gz -2 ", args$illumina, "_R2_001.fastq.gz | samtools view -bS - | samtools sort -o ", pilondir, "/btbam/", args$canupre, ".sorted.bam"))
    system(paste0("samtools index ", pilondir, "/btbam/", args$canupre, ".sorted.bam"))

    
    ##Use pilon to fix nanopore contigs with alginments
    system(paste0("java -Xmx30G -jar ~/software/pilon/pilon-1.22.jar --threads 12 --changes --tracks --genome ", genome, " --frags ", pilondir, "/btbam/", args$canupre, ".sorted.bam", " --outdir ", pilondir, " --output ", args$canupre))
    system(paste0('cp ', pilondir,'/', args$canupre, '.fasta ', pilondir, '/btidx/')) 
    genome=file.path(pilondir, 'btidx', paste0(args$canupre,".fasta"))
    system(paste0("sed -i -e 's/_pilon//g' ", genome))

    print(paste0('Finished round ', as.character(i)))
    print('=================================================================================================================================================================')
    
}

    
