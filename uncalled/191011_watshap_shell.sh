#!/bin/bash 

main=/kyber/Data/Nanopore/Analysis/gilfunk/190824_allGuidez_take2_GM12878
samp=180924_allGuideMin

### this assumes you already have  a vcf file 
# this script then takes those variants in vcf and makes a phased vcf ( '|'   instead of '\' ) 
 # i dont know how it does so if the variants are not linked .. lets see what the output looks like 
# ahh, i see now.. 
#it seems it uses blocks and reads in those blocks to do phasing ..  and then assigns releative to those other variants in the blocK !! 

in_vcf=$main/varIANT_calls/nanoP_190824_allGuidez_take2_GM12878_VARmin20/*commonfromONTARGvariants.vcf
in_bam=$main/prim_onTarg/190824_allguides_take2_gm12878_onTarg_prim.bam

out_dir=$main/watshap_nanopolish
#$main/varTISSUEnp/$samp
phased_vcf_out=$out_dir/${samp}_npSNPs_phased.vcf
tagged_bam=$out_dir/${samp}_tagged_WH.bam

ref=/mithril/Data/NGS/Reference/human38/GRCH38.fa

  # .. you need the 'ignore read groups'  or else it will look for sample tags in the bam 
   # also rn i did not provide a reference,  but you could --  unclear to me how that would affec the output 
 # apparently i did include a reference now tho 

# # THIS SECTION is only needed if your vcf is NOT phased --  
  # #  
if false ; then
  echo 'just you wait your pretty little self there while watsHap does its thing'
  echo 'making phased vcf'
  /home/gilfunk/.local/bin/whatshap phase  --full-genotyping  --indels --reference $ref --ignore-read-groups -o $phased_vcf_out $in_vcf $in_bam 

  echo 'compressing vcf with bgzip' 
  bgzip ${phased_vcf_out}
  echo 'indexing vcf with tabix'
  tabix -p vcf ${phased_vcf_out}.gz

fi 

### this next bit takes a 'phased' vcf, then looks at the reads once again and assigns them to a read group 

if false  ; then 
  echo 'making a new bam with haplotype tags'
  /home/gilfunk/.local/bin/whatshap haplotag  --ignore-read-groups  -o  $tagged_bam  ${phased_vcf_out}.gz  $in_bam
  echo 'indexing tagged bam file'
   samtools index $tagged_bam

   echo 'starting splitting of that bam bieebbyyyy!'  
   python helper_scripts/v2_split_BAM_by_tag.py  -b $tagged_bam  -o $out_dir  -n $samp

   echo 'indexing split  bam fileSSS'
   samtools index $out_dir/*hap1.bam
   samtools index $out_dir/*hap2.bam
fi


## this last bit is for getting fastqs for each allele for assembly and consensus 
if false ; then
   echo 'converting BAMs to Fastq'
   samtools fastq $out_dir/*hap1.bam > $out_dir/${samp}_hap1.fastq
   samtools fastq $out_dir/*hap2.bam > $out_dir/${samp}_hap2.fastq
fi


