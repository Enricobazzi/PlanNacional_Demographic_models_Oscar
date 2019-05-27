#!/bin/bash
#SBATCH -t 0
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH ​--mail-type=begin
#SBATCH ​--mail-type=end
#SBATCH ​--mail-user=enricobazzical@gmail.com


###########################################
## GATK 4.1.1.0 HaplotypeCaller launcher ##
###########################################

# This program runs GATK 4.1.1.0 HaplotypeCaller with standard settings for per-sample
# GVCF generation on the finis terrae II server of CESGA.

# The sample name must be defined while launching the script as such:

# ./cesga_Haplotypecaller.sh <sample>

module load gatk/4.1.1.0

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/test/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# Input File:
INbam=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_bams/"$1"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam

# Output Files:
OUTgvcf=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_gvcfs/"$1"_cat_ref.g.vcf.gz

##################################
## GATK 4.1.0.0 HaplotypeCaller ##
##################################

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I $INbam \
   -O $OUTgvcf \
   --native-pair-hmm-threads 24 \
   -ERC GVCF
