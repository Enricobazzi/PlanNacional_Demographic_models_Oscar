#!/bin/bash

###########################################
## GATK 4.1.0.0 HaplotypeCaller launcher ##
###########################################

# This program runs GATK 4.1.0.0 HaplotypeCaller with standard settings for per-sample
# GVCF generation for one specific chromosome.

# The sample name and chromosome name must be defined while launching the script as such:

# ./sample_chr_Haplotypecaller.sh <sample> <chromosome>

# A specific BED with the chromosome name must be previously generated in desired location
# with specific name (see Parallel_Calling_pipeline.md file for further information on their generation)

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# Input File:
INbam=/home/ebazzicalupo/CatRef_bams/"$1"_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam

# chromosome BED file:
BED=/home/ebazzicalupo/CatGenome_CHR_BEDs/"$2"_CHR_coordinates.bed

# Output Files:
OUTgvcf=/home/ebazzicalupo/CatRef_gvcfs/"$1"_"$2"_cat_ref.g.vcf.gz

##################################
## GATK 4.1.0.0 HaplotypeCaller ##
##################################

/opt/gatk-4.1.0.0/gatk HaplotypeCaller  \
   -R $REF \
   -I $INbam \
   -O $OUTgvcf \
   -L $BED \
   -ERC GVCF
