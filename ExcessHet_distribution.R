library(tidyverse)

Excess_het <- read_table("/Users/enricobazzicalupo/Documents/PlanNacional/gatk_stats_distributions/ExcessHet.table",
                   col_names = F) %>%
  rename("ExcessHet" =  X1) %>% as.data.frame()

freq_table_DF <- as.data.frame(table(as.numeric(Excess_het$ExcessHet)))

ggplot(Excess_het, aes(x=ExcessHet)) +
  geom_histogram(bins = 30) +
  scale_x_continuous(breaks = 0:60*5, limits = c(0, 60)) +
  scale_y_continuous(limits = c(0, 100000))




  