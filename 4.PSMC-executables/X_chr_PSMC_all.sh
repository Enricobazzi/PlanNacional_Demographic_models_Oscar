#!/bin/bash

## X chromosome pseudodiploid PSMC analysis for estimating divergence times ##

cd ~/PSMC

# All pair combinations are too many
# One per spp (two for ll -> east + west) will be better to start with:

IND_LIST=($(ls consensus_fq/*.fq | rev | cut -d'/' -f1 | rev | grep -v "X_chr" | cut -d'.' -f1 | grep -v "-" | grep -E "c_lc_zz_0003|c_ll_ki_0090|c_ll_ya_0146|c_lp_sm_0138|c_lr_fl_0005"))

## 1. Get fq of X chromosome only
# I use seqtk subseq:

for ind in ${IND_LIST[@]}
 do

  echo "Extracting X_chr.fq of ${ind}"
  seqtk subseq consensus_fq/${ind}.fq ~/CatGenome_CHR_BEDs/X_CHR_coordinates.bed \
  > consensus_fq/${ind}.X_chr.fq

done

# Then work on Pairs of individuals

START=0 # arrays start counting from 0
END=($(echo ${IND_LIST[@]} | wc -w)) # after put $END -1 because arrays start counting from 0

for i in $(seq $(($START)) $(($(($END))-1)));
 do
  for j in $(seq $(($i + 1)) $(($(($END))-1)));
   do

   ## 2. Merge X_chr.fq of the two individuals
   # I use seqtk mergefa
   echo "Merging ${IND_LIST[${i}]} and ${IND_LIST[${j}]}"
   seqtk mergefa consensus_fq/${IND_LIST[${i}]}.X_chr.fq consensus_fq/${IND_LIST[${j}]}.X_chr.fq \
   > consensus_fq/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.fq

   ## 3. Convert my pseudohaploid sample fq to PSMC input format
   # I use fq2psmcfa
   echo "converting ${IND_LIST[${i}]}-${IND_LIST[${j}]} to PSMC format"
   /opt/psmc/utils/fq2psmcfa consensus_fq/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.fq \
   > inputfiles_psmcfa/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.psmcfa

   ## 4. Run PSMC on the pseudohaploid
   echo	"Running PSMC on ${IND_LIST[${i}]}-${IND_LIST[${j}]}"
   psmc -p "4+25*2+4+6" -o output_psmc/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.psmc \
   inputfiles_psmcfa/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.psmcfa

   ## 5. Paint results
   echo	"Painting ${IND_LIST[${i}]}-${IND_LIST[${j}]} results"
   /opt/psmc/utils/psmc_plot.pl -u 1.6e-08 -g 5 -Y 20 plots/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.plot \
   output_psmc/${IND_LIST[${i}]}-${IND_LIST[${j}]}.X_chr.psmc

  done
done
