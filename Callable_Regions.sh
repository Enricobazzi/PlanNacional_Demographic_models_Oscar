#!/bin/bash

# We decided to filter our bams instead of our VCF for depth and missingness (indirectly)
# I will use bedtools genomecov to create a bed of the depth at each window in the genome.
# This way I can filter this bed for the desired values of depth (>5X and <MAXdepth) and create
# a BED of the "callable regions of the genome"

cd ~/8samples_callable

# To run it in parallel I developed a small script (callregions.sh) that will be run in different screens 
# at the same time on the genomics-b server:

#!/bin/bash
samplepath=($(echo "${1}"))
sample=($(echo "${samplepath}" | rev | cut -d "/" -f1 | rev | cut -d'_' -f1-4))

bedtools genomecov \
-ibam ${samplepath}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam -bga \
> ${sample}.coverage.bed

# To run the script I just need to give it the sample path, in a loop:
samplearray=(contemporary_data/c_ll_vl_0112 contemporary_data/c_ll_ya_0146 contemporary_data/c_lp_sm_0138 contemporary_data/c_lp_sm_0140 c_lc_zz_0001 c_lc_zz_0003 c_lr_zz_0001 c_lr_nm_0006)

for sample in ${samplearray[@]}
 do 
  samplename=($(echo "${sample}" | rev | cut -d "/" -f1 | rev | cut -d'_' -f1-4))
  echo "${samplename}"
  screen -dmS callable_${samplename}  sh -c "/home/ebazzicalupo/callregions.sh /GRUPOS/grupolince/CatRef_bams/${sample}; exec /bin/bash"
done

# Recalculate the limits based on individual distributions of the samples I want to use:

####
# IN R:
library(tidyverse)
# I want to calculate upper limits of depth for the samples we decided on using with Oscar in order 
# create a callable regions bed file for each sample.

# For LP:
# I see by comparing the depth files that the 4th column is c_lp_sm_0138 and the 5th is c_lp_sm_0140

# For LL:
# I see by comparing the depth files that the 3rd column is c_ll_vl_0112 and the 4th is c_ll_ya_0146

# Basically they go in the order they are in the bamlist (for future reference)
# So I can re-analyze the distributions from the depth files I already have on my laptop

wd_input <- "/Users/enricobazzicalupo/Documents/allsamples_depths/"

# Create a list of the sample files, sample names and the corresponding column
sample_files <- c("LCmacro", "LCMurphy", "LLmacro", "LLmacro", "LPpgenoma", "LPpgenoma", 
                  "LRJan", "LRmacro")

samples <- c("c_lc_zz_0001", "c_lc_zz_0003", "c_ll_vl_0112", "c_ll_ya_0146",
             "c_lp_sm_0138", "c_lp_sm_0140", "c_lr_nm_0006", "c_lr_zz_0001")

columns <- c(3, 3, 5, 6, 6, 7, 3, 3)

# Create an Empty dataframe to save values for each dataset for the final table
depth_per_sample <- data.frame()

################################################
##### Depth Calculations, Graphs and Table #####
################################################
i=1
# For every Sample file:
for (i in 1:length(samples)){
  
  sample_file <- sample_files[i]
    
  sample <- samples[i]
   
  column <- columns[i]
  
  # Import the sample file table
  input.depth <- read_delim(paste0(wd_input,sample_file, ".200x100kbp.masked.depth"), 
                            col_names = F, delim = '\t')
  
  # create sample's depth table
  c1 <- colnames(input.depth)[column]
  sample_depth <- data.frame(input.depth[,column])
  names(sample_depth)[names(sample_depth) == c1] <- "depth"
  
  # Create a frequency table of the values of the Total column
  freq_table_DF <- as.data.frame(table(sample_depth$depth))
  
  # Define the functions for mean and standard deviation, and define maximum and minimum depth based on different criteria
  mean_folds = 0.95
  
  my_mean_DF <- mean(sample_depth$depth)
  my_sd_DF <- sd(sample_depth$depth)
  
  maxDepth_DF_meanfolds = my_mean_DF + (mean_folds * my_mean_DF)
  minDepth_DF_meanfolds  = my_mean_DF - (mean_folds * my_mean_DF)
  
  maxDepth_DF_SD = my_mean_DF + (my_sd_DF * 1.5)
  minDepth_DF_SD = my_mean_DF - (my_sd_DF * 1.5)
  
  maxDepth_MEAN <- my_mean_DF * 2
  minDepth_5x <- 5
  
  # Calculate percentage of variants left after filtering for min Depth = 5x n. individuals, and max Depth = mean + 1.5 sd
  PercentLeft <- length(which(sample_depth$depth > minDepth_5x & sample_depth$depth < maxDepth_DF_SD)) / nrow(sample_depth) *100
  
  # Add dataset information to Dataframe
  depth_per_sample <- rbind(depth_per_sample,
                            data.frame(Sample = sample,
                                       mean = my_mean_DF, sd = my_sd_DF, prcentPosLeft = PercentLeft,
                                       maxDepthMF = maxDepth_DF_meanfolds, minDepthMF = minDepth_DF_meanfolds,
                                       maxDepthSD = maxDepth_DF_SD, minDepthSD = minDepth_DF_SD,
                                       maxDepthMEAN = maxDepth_MEAN, minDepth5x = minDepth_5x))
  
  # Draw and save a graph of the distribution of depth values, with upper and lower depth limits
  ggplot(freq_table_DF, aes(x = as.numeric(Var1), y = Freq)) +
    geom_bar(stat = "identity", color = "black") +
    scale_x_continuous(breaks = 0:250*10, limits = c(0, maxDepth_DF_SD*1.5)) +
    # scale_x_discrete(limits = c(0, maxDepth_DF*1.5)) +
    scale_y_continuous(expand=c(0,0)) +
    geom_vline(xintercept=maxDepth_DF_SD,linetype="dashed", size=0.5) +
    #geom_vline(xintercept=minDepth_DF_SD,linetype="dashed", size=0.5) +
    #geom_vline(xintercept=minDepth_DF_meanfolds, colour ="grey", linetype="dashed", size=0.5) +
    #geom_vline(xintercept=maxDepth_DF_meanfolds, colour ="grey", linetype="dashed", size=0.5) +
    #geom_vline(xintercept=minDepth_5x, colour ="red", linetype="dashed", size=0.5) +
    #geom_vline(xintercept=maxDepth_MEAN, colour ="red", linetype="dashed", size=0.5) +
    theme_classic() +
    theme(text = element_text(size=10))
  ggsave (filename = (paste0("graph_",sample,".depth.pdf")), path = wd_input)
}

# Print the table to a file
write.table(x = depth_per_sample,file = paste0(wd_input,"8samples_depth.csv"),quote=FALSE, col.names = T, row.names = FALSE, sep= ",")

####

# Upload the table to genomics-b:
scp /Users/enricobazzicalupo/Documents/allsamples_depths/8samples_depth.csv ebazzicalupo@genomics-b.ebd.csic.es:~

# Now to assign "not-callable" status to <5X and >maxdepth regions

for sample in $(ls *.coverage.bed | cut -d'.' -f1)
 do 
  echo ${sample}
  max=($(grep "${sample}" 8samples_depth.csv | cut -d',' -f7))
  awk -v max="${max}" '{FS="\t"; OFS="\t"; if ($4 < 6 || $4 > max) print $1, $2, $3, "not-callable"; else print $1, $2, $3, "callable";}' ${sample}.coverage.bed \
  > ${sample}.coverage.defined.bed
done

# I can now remove the not-callable regions and also the repetitive/low-mappability ones,
# to obtain the callable regions of each sample:

for sample in $(ls *.coverage.defined.bed | cut -d'.' -f1)
 do 
  echo ${sample}
  grep -v "not-callable" ${sample}.coverage.defined.bed | 
  bedtools subtract -a - -b /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Masked_Regions.bed |
  bedtools merge -i - \
  > ${sample}.callable.bed
done

for sample in $(ls *.coverage.defined.bed | cut -d'.' -f1)
 do 
  bedtools merge -i ${sample}.callable.bed > ${sample}.callable.collapsed.bed
done

# To check length of sequence in BED
# awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' c_lc_zz_0001.callable.collapsed.bed

# Intersect all of the BED files
multiIntersectBed -i c_lc_zz_0001.callable.collapsed.bed \
c_lc_zz_0003.callable.collapsed.bed \
c_ll_vl_0112.callable.collapsed.bed \
c_ll_ya_0146.callable.collapsed.bed \
c_lp_sm_0138.callable.collapsed.bed \
c_lp_sm_0140.callable.collapsed.bed \
c_lr_nm_0006.callable.collapsed.bed \
c_lr_zz_0001.callable.collapsed.bed \
> 8samples.callable.collapsed.intersect.bed

# Only regions overlapping and callable in all samples
grep "1,2,3,4,5,6,7,8" 8samples.callable.collapsed.intersect.bed > 8samples.callable.collapsed.intersect.all.bed

# Get VCF of 8samples from allsamples VCF
samplesARRAY=($(ls *.coverage.defined.bed | cut -d'.' -f1))

/opt/gatk-4.1.0.0/gatk SelectVariants \
  -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  -V ~/CatRef_vcfs/allsamples_cat_ref.filter5.vcf \
  $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
  -O 8samples_cat_ref.filter5.vcf


# Now we can filter the VCF 
bedtools subtract -a 8samples_cat_ref.filter5.vcf \
  -b 8samples.callable.collapsed.intersect.all.bed \
  -header > 8samples_cat_ref.filter-callable.vcf
  
# Copy to laptop to share with Oscar
cd Documents/Oscar_v3_vcfs
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/8samples_callable/8samples.callable.collapsed.intersect.all.bed.gz .
scp ebazzicalupo@genomics-b.ebd.csic.es:/home/ebazzicalupo/8samples_callable/8samples_cat_ref.filter-callable.vcf.gz .
