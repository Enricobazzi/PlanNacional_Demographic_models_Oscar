#!/bin/bash
#SBATCH -t 2-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com


###########################################
## GATK 4.1.1.0 HaplotypeCaller launcher ##
###########################################

# This program runs GATK 4.1.1.0 HaplotypeCaller with standard settings for per-sample per-chromosome
# GVCF generation on the finis terrae II server of CESGA.

# The sample name and chromosome must be defined while launching the script as such:

# ./cesga_Haplotypecaller.sh <sample> <chromosome>

module load gatk/4.1.1.0

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/test/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# Input File:
INbam=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_bams/"$1"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam

# Output Files:
OUTgvcf=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_gvcfs/"$1"_"$2"_cat_ref.g.vcf.gz

# chromosome BED file:
BED=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatGenome_CHR_BEDs/"$2"_CHR_coordinates.bed


##################################
## GATK 4.1.1.0 HaplotypeCaller ##
##################################

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I $INbam \
   -O $OUTgvcf \
   -L $BED \
   --native-pair-hmm-threads 24 \
   -ERC GVCF
