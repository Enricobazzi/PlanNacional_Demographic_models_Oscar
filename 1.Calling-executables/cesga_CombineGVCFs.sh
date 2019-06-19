#!/bin/bash
#SBATCH -t 2-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

####################################################
START=$(date)
echo "CombineGVCFs SCRIPT for $1 starting : $START"
####################################################

########################################
## GATK 4.1.1.0 CombineGVCFs launcher ##
########################################

# This program runs GATK 4.1.1.0 CombineGVCFs with standard settings for per-chromosome
# GVCF generated on the finis terrae II server of CESGA.

# The sample name must be defined while launching the script as such:

# ./cesga_CombineGVCFs.sh <chromosome>

module load gatk/4.1.1.0

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/test/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# Input File:
INgvcfsARRAY=($(ls /mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_gvcfs/*_"$1"_cat_ref.g.vcf.gz))

# Output Files:
OUTgvcf=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatRef_gvcfs/"$1"_cat_ref_joint.g.vcf.gz

# chromosome BED file:
BED=/mnt/lustre/scratch/home/csic/ebd/jg2/test/CatGenome_CHR_BEDs/"$1"_CHR_coordinates.bed


###############################
## GATK 4.1.1.0 CombineGVCFs ##
###############################

gatk --java-options "-Xmx4g" CombineGVCFs \
   -R $REF \
   $(for i in ${INgvcfsARRAY[@]}; do echo "--variant ${i}";done) \
   -L $BED \
   -O $OUTgvcf

####################################################
END=$(date)
echo "CombineGVCFs SCRIPT for $1 ended : $END"
####################################################
