#!/bin/bash
#SBATCH -t 9:00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

#############################################################
START=$(date)
echo "Missingness_Filter_6 SCRIPT for $1 starting : $START"
#############################################################

#############################
## Samtools Depth launcher ##
#############################

# With this script I want to calculate depth at each position of various bamlist files
# using samtools depth. Depth at all positions will be calculated (-a) within the
# regions randomly selected before (-b) (see 2.Variant_Filtering.md for more detail).
# This will be run in a loop for all bamlists.

# Samtools version 1.9 will be used.

module load samtools

###################################
## VARIABLE and PATHS definition ##
###################################

# Array of bamlist files
BAMlistArray=($(ls $LUSTRE/test/CatRef_bams/*.bamlist | rev | cut -d'/' -f 1 | rev))

#####################################
## Calculating Depth with SAMtools ##
#####################################

# Loop of Samtools depth calculations for each bamlist
for i in ${BAMlistArray[@]}
  do
  echo "Calculating depth for $i"
  samtools depth -a -b $LUSTRE/test/FilterTrials/Felis_catus.200x100kbp.masked.genome.bed \
  -f $LUSTRE/test/CatRef_bams/"$i" \
  > $LUSTRE/test/FilterTrials/"$i".200x100kbp.masked.depth
done

###########################################################
END=$(date)
echo "Missingness_Filter_6 SCRIPT for $1 ended : $END"
###########################################################
