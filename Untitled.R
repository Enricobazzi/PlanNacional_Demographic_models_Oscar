#############
## TREEMIX ##
#############

source("~/plotting_funcs.R")
plot_tree("~/try1/treemix_out")
plot_tree("~/try2/treemix_out")

poporder <- c("lc", "ll", "lp", "lr")
plot_resid("~/try1/treemix_out", "~/try1/poporder")
plot_resid("~/try2/treemix_out", "~/try1/poporder")
