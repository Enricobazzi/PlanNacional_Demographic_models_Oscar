#!/bin/bash
#SBATCH -t 2-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

#############################################################
START=$(date)
echo "Depth_Filter_7 SCRIPT for $1 starting : $START"
#############################################################

##############################
## Depth Filtering launcher ##
##############################

# With this script I want to apply a filter based on depth to my different
# datasets. As explained in 2.Variant_Filtering.md my datasets are:

# Lynx lynx - one dataset :
#   group 1 - cr_0212, ki_0090, vl_0112, ya_0146 with MACROGEN : LLmacro

# Lynx rufus - three datasets :
#   group 1 - c_lr_fl_0005 from Murphy : LRMurphy
#   group 2 - c_lr_nm_0006 from Janecka : LRJan
#   group 3 - c_lr_zz_0001 with MACROGEN : LRmacro

# Lynx canadiensis - two datasets :
#   group 1 - c_lc_zz_0001 with MACROGEN : LCmacro
#   group 2 - c_lc_zz_0003 from Murphy : LCMurphy

# Lynx pardinus - two datasets :
#   group 1 - c_lp_do_0153, c_lp_do_0173, c_lp_do_0443, c_lp_sm_0138, c_lp_sm_0140, c_lp_sm_0185, c_lp_sm_0186,
#              c_lp_sm_0298, c_lp_sm_0359, h_lp_do_0007 from Proyecto Genoma : LPpgenoma
#   group 2 - c_lp_sm_0221 from CANDILES : LPcandiles

# After dividing my global VCF into single dataset VCFs (with only individuals of the same dataset),
# I will apply a filter to each dataset based on the maximum and minimum values
# of depth calculated with the R script depth_loop.R. These values are saved into a table called
# depth_per_sample.csv stored in FilterTrials directory.
# I will finally filter out from the global VCF the variants filtered out in the single
# dataset VCFs.

# This script will use the following softwares:

# BEDtools 2.28.0
# BCFtools 1.9
# GATK 4.1.0.0

module load bedtools
module load bcftools
module load gatk/4.1.1.0

# # The VCF file name (without the extensions) must be defined while launching
# the script as such:

# ./Depth_Filter_7.sh <VCFfilename>

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/mnt/lustre/scratch/home/csic/ebd/jg2/test/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa

# AllIndividuals input VCF
INVCF=$LUSTRE/test/CatRef_vcfs/"$1".filter6.vcf

# VCF Directory
OUTdir=$LUSTRE/test/CatRef_vcfs

# Depth per sample Table generated in R (depth_loop.R)
DPStable=$LUSTRE/test/FilterTrials/depth_per_sample.csv

# Create a copy of the VCF file with a new name that will be used to filter out
# variants with excessively low/high depth:
cp $INVCF $OUTdir/"$1".filter7.vcf
OUTVCF=$OUTdir/"$1".filter7.vcf

# Create a log file for keeping track of variant numbers in the various steps
echo "Per dataset Depth variant filtering:" > $LUSTRE/test/depth.filtering.variants.log
depthLOG=$LUSTRE/test/depth.filtering.variants.log

#############################
## Applying filters - LOOP ##
#############################

# Report starting number of variants
startVARs=($(grep -v "#" $INVCF | wc -l))
echo "starting variants : $startVARs" >> $depthLOG

# List of all datasets in an array:
datasetARRAY=($(ls $LUSTRE/test/CatRef_bams/*.bamlist | rev | cut -d'/' -f1 | rev | cut -d'.' -f1))

# For each dataset:
for i in ${datasetARRAY[@]}
  do

    # (1) Extract VCF of individuals of dataset
    samplesARRAY=($(cat $LUSTRE/test/CatRef_bams/${i}.bamlist | rev | cut -d'/' -f1 | rev | cut -d'_' -f1-4))

    gatk SelectVariants \
    -R $REF \
    -V $INVCF \
    $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
    -O $OUTdir/${i}.vcf

    # Minimum and Maximum depth values (in this case a min of 5x per sample = column 11, and a max of mean+1.5*sd = column 8)
    # see depth_loop.R for more detail on each column value
    max=$(grep ${i} $DPStable | cut -d',' -f8)
    min=$(grep ${i} $DPStable | cut -d',' -f11)
    echo "Maximum depth of ${i} is ${max}, Minimum depth is ${min}"

    # (2) extract the excessive missingness variant with BCFtools filter
    echo "extracting excessively low/high depth variants from $i VCF"
    bcftools filter -i "INFO/DP < ${min} || INFO/DP > ${max}" -Ov $OUTdir/${i}.vcf \
    > $OUTdir/${i}.applydepthfilter.vcf

    # (3) filter the excessively missing variants from the new VCF file of all samples
    echo "subtracting excessively low/high depth variants of $i from output VCF"
    bedtools subtract -a $INVCF \
    -b $OUTdir/${i}.applydepthfilter.vcf -header \
    > tmp && mv tmp $OUTVCF

    # Report number of variants in dataset
    datasetVARs=($(grep -v "#" $OUTdir/${i}.vcf | wc -l))
    echo "${i} variants : $datasetVARs" >> $depthLOG
    # Report number of variants filtered out
    filteredVARs=($(grep -v "#" $OUTdir/${i}.vcf | wc -l))
    echo "number of variants filtered from ${i} : $filteredVARs" >> $depthLOG

done

# Report final number of variants
finalVARs=($(grep -v "#" $OUTVCF | wc -l))
echo "final variants : $finalVARs" >> $depthLOG


###########################################################
END=$(date)
echo "Depth_Filter_7 SCRIPT for $1 ended : $END"
###########################################################
