#!/bin/bash

############################
## Parallel Caller Script ##
############################

# THIS VERSION IS FOR 2 SAMPLES ONLY : lp_sm_0359 & ll_cr_0212 (see array)

# This program will create a screen and execute the sample_chr_Haplotypecaller.sh
# script for each sample in an array of sample BAMs, for a single chromosome defined by
# an element of the array of chromosome BEDs. The BED array element number must be defined
# while launching the script as such:

# ./chr_parallel_caller.sh <N>

# where N is the BED array element number. Here is the list of all elements in the BED array
# with their corresponding number:

# 0 = A1
# 1 = A2
# 2 = A3
# 3 = B1
# 4 = B2
# 5 = B3
# 6 = B4
# 7 = C1
# 8 = C2
# 9 = D1
# 10 = D2
# 11 = D3
# 12 = D4
# 13 = E1
# 14 = E2
# 15 = E3
# 16 = F1
# 17 = F2
# 18 = MT
# 19 = rest
# 20 = X

# More information on the BEDs and the script can be found in the 1.Parallel_Calling_pipeline.md file

# CAREFUL: sample_chr_Haplotypecaller.sh MUST be in the same folder where this is launched!!!


###################################
## VARIABLE and PATHS definition ##
###################################

# Array of BED files
BEDarray=($(ls /home/ebazzicalupo/CatGenome_CHR_BEDs/*.bed | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1 | uniq))

# Array of BAM files (samples)
SAMPLEarray=($(ls /home/ebazzicalupo/CatRef_bams/*.bam | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1,2,3,4 | uniq | grep -E "lp_sm_0359|ll_cr_0212"))

##########################################
## Creation of a Screen for each sample ##
##########################################

for bed in ${BEDarray[$1]}
  do
    echo ${bed}
    for sample in ${SAMPLEarray[@]}
      do
        screen -dmS "${sample}_${bed}_calling"
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "script ${sample}_${bed}_calling.log\n"
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "time ./try_sample_chr_Haplotypecaller.sh ${sample} ${bed}; exec bash\n"
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "exit\n"
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "exit\n"
    done
done
