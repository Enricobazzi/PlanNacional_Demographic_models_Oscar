#!/bin/bash

############################
## Parallel Caller Script ##
############################

# This program will create a screen and execute the sample_chr_Haplotypecaller.sh
# script for each sample in an array of samples, for a single chromosome defined by
# an element of the array of chromosome BEDs. The BED array element number must be defined
# while launching the script as such:

# ./chr_parallel_caller.sh <N>

# where N is the BED array element number. Here is the list of all elements in the array
# with their corresponding number:

# 1 = A1
# 2 = A2
# 3 = A3
# 4 = B1
# 5 = B2
# 6 = B3
# 7 = B4
# 8 = C1
# 9 = C2
# 10 = D1
# 11 = D2
# 12 = D3
# 13 = D4
# 14 = E1
# 15 = E2
# 16 = E3
# 17 = F1
# 18 = F2
# 19 = MT
# 20 = rest
# 21 = X

# More information on the BEDs and the script can be found in the 1.Parallel_Calling_pipeline.md file


###################################
## VARIABLE and PATHS definition ##
###################################

# Array of BED files
BEDarray=($(ls /home/ebazzicalupo/CatGenome_CHR_BEDs/*.bed | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1 | uniq))

# Array of BAM files (samples)
SAMPLEarray=($(ls /home/ebazzicalupo/CatRef_bams/*.bam | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1,2,3,4 | uniq))

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
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "time ./sample_chr_Haplotypecaller.sh ${sample} ${bed}; exec bash\n"
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "exit\n"
        screen -S "${sample}_${bed}_calling" -p 0 -X stuff "exit\n"
    done
done
