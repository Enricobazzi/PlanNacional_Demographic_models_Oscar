##############################
##### DEPTH CALCULATIONS #####
##############################

# This R script will be used to analyze the depth per position tables of the different datasets
# generated with SAMTOOLS depth.

# The different tables contain the information of of the depth at each position for each individual,
# for a randomly generated subset of the genome. After summing the depth of all the individuals of the
# dataset, we can make a frequency table of depth values for the dataset.

# Defining different upper and lower limits to the distribution have been tried:
# based on mean (upper and lower 95% of mean, or a maximum 2 times the mean);
# standard deviation (upper and lower 1.5 times the standard deviation);
# and fixed parameters (a minimum of 5x depth per sample).

# The graphs generated help visualize where these limits fall within the distribution and decide which
# ones better fit our dataset.

# A table with the different values will be generated and printed into a file (),
# so this file can be later used in Bash to retract the different values for filtering the dataset.

# The final criteria chosen for the limits was ...

##########################
##### Load Libraries #####
##########################

library(readr)
library(dplyr)
library(ggplot2)

######################################
##### Import and Prepare Dataset #####
######################################

# Define the input Directory
wd_input <- "/Users/enricobazzicalupo/Documents/PlanNacional/gatk_stats_distributions/"

# Create a list of the sample files' names (the SAMTOOLS depth output files)
sample_files <- list.files(wd_input, pattern="*.depth$")

# Create an Empty dataframe to save values for each dataset for the final table
depth_per_sample <- data.frame()

################################################
##### Depth Calculations, Graphs and Table #####
################################################

# For every Sample file:
for (i in 1:length(sample_files)){

  # Import the sample file table
  input.depth <- read_delim(paste0(wd_input,sample_files[[i]]), col_names = F, delim = '\t')

  # Add a column (Total) which is the sum of all the depth columns (from the third to the last)
  input.depth$Total <- rowSums(input.depth[,3:ncol(input.depth)])

  # Create a frequency table of the values of the Total column
  freq_table_DF <- as.data.frame(table(input.depth$Total))

  # Define the functions for mean and standard deviation, and define maximum and minimum depth based on different criteria
  mean_folds = 0.95

  my_mean_DF <- mean(input.depth$Total)
  my_sd_DF <- sd(input.depth$Total)

  maxDepth_DF_meanfolds = my_mean_DF + (mean_folds * my_mean_DF)
  minDepth_DF_meanfolds  = my_mean_DF - (mean_folds * my_mean_DF)

  maxDepth_DF_SD = my_mean_DF + (my_sd_DF * 1.5)
  minDepth_DF_SD = my_mean_DF - (my_sd_DF * 1.5)

  maxDepth_MEAN <- my_mean_DF * 2
  minDepth_5x <- 5 * (ncol(input.depth) - 3)

  # Calculate percentage of variants left after filtering for min Depth = 5x n. individuals, and max Depth = mean + 1.5 sd
  PercentLeft <- length(which(input.depth$Total > minDepth_5x & input.depth$Total < maxDepth_DF_SD)) / nrow(input.depth) *100

  # Define Dataset Name:
  population=unlist(strsplit(basename(sample_files[[i]]),"[.]"))[1]

  basename(sample_files[[i]])
  # Add dataset information to Dataframe
  depth_per_sample <- rbind(depth_per_sample,
                            data.frame(pop = population, nindividuals = (ncol(input.depth) - 3),
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
 ggsave (filename = (paste0("graph_",basename(sample_files[[i]]),".pdf")), path = wd_input)
}

# Print the table to a file
write.table(x = depth_per_sample,file = paste0(wd_input,"depth_per_sample.csv"),quote=FALSE, col.names = T, row.names = FALSE, sep= ",")
