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
# depth_per_sample.csv stored in XXX directory.
# I will finally filter out from the global VCF the variants filtered out in the single
# dataset VCFs.

# This script will use the following softwares:

# BEDtools 2.28.0
# BCFtools 1.9

module load bedtools
module load bcftools

# # The VCF file name (without the extensions) must be defined while launching
# the script as such:

# ./Depth_Filter_7.sh <VCFfilename>

#####################################
## Applying filters - Preparations ##
#####################################

# List all sample names in a .namelist file:
ls $LUSTRE/test/CatRef_bams/*.bam | rev | cut -d'/' -f1 | rev | cut -d '_' -f1-4 | sort -u \
> $LUSTRE/test/CatRef_bams/all-samples.namelist

# Declare dataset elements
declare -A BARCODEID=(["c_lc_zz_0001"]="LCmacro" ["c_lc_zz_0003"]="LCMurphy" ["c_ll_cr_0212"]="LLmacro" \
["c_ll_ki_0090"]="LLmacro" ["c_ll_vl_0112"]="LLmacro" ["c_ll_ya_0146"]="LLmacro" ["c_lp_do_0153"]="LPpgenoma" \
["c_lp_do_0173"]="LPpgenoma" ["c_lp_do_0443"]="LPpgenoma" ["c_lp_sm_0138"]="LPpgenoma" ["c_lp_sm_0140"]="LPpgenoma" \
["c_lp_sm_0185"]="LPpgenoma" ["c_lp_sm_0186"]="LPpgenoma" ["c_lp_sm_0298"]="LPpgenoma" ["c_lp_sm_0359"]="LPpgenoma" \
["h_lp_do_0007"]="LPpgenoma" ["c_lr_fl_0005"]="LRMurphy" ["c_lr_nm_0006"]="LRJan" ["c_lr_zz_0001"]="LRmacro" \
["c_lp_sm_0221"]="LPcandiles")

# List individuals in an array (for loop):
individualsARRAY=($(cat $LUSTRE/test/CatRef_bams/all-samples.namelist))

for i in ${individualsARRAY[@]}
  do
    echo ${BARCODEID["${i}"]}
done



###########################################################
END=$(date)
echo "Depth_Filter_7 SCRIPT for $1 ended : $END"
###########################################################
