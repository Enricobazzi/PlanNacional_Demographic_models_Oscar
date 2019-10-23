#############
## TREEMIX ##
#############

source("~/plotting_funcs.R")
plot_tree("~/try1/treemix_out")
plot_tree("~/try2/treemix_out")

plot_resid("~/try1/treemix_out", "~/try1/poporder")
plot_resid("~/try2/treemix_out", "~/try1/poporder")

# Treemix Without pre-defined topology and no Outgroup
for (i in 1:10){
  plot_tree(paste0("~/Documents/Treemix_results/m_",i,"/treemix_out"))
  plot_resid(paste0("~/Documents/Treemix_results/m_",i,"/treemix_out"), "~/Documents/Treemix_results/poporder")
  
}

# Treemix Without pre-defined topology and Cat as Outgroup
for (i in 1:10){
  plot_tree(paste0("~/Documents/Treemix_outgroup_results/m_",i,"/treemix_out"))
  plot_resid(paste0("~/Documents/Treemix_outgroup_results/m_",i,"/treemix_out"), "~/Documents/Treemix_results/poporder_cat")
  
}

# Treemix WITH pre-defined topology and LR Outgroup
for (i in 1:10){
  plot_tree(paste0("~/Documents/Treemix_results/m_topo_",i,"/treemix_out"))
  plot_resid(paste0("~/Documents/Treemix_results/m_topo_",i,"/treemix_out"), "~/Documents/Treemix_results/poporder")
  
}

# Treemix WITH pre-defined topology and Cat as Outgroup
for (i in 1:10){
  plot_tree(paste0("~/Documents/Treemix_outgroup_results/m_topo_",i,"/treemix_out"))
  plot_resid(paste0("~/Documents/Treemix_outgroup_results/m_topo_",i,"/treemix_out"), "~/Documents/Treemix_results/poporder_cat")
  
}

# Treemix WITH pre-defined topology and Cat as Outgroup K=50
for (i in 1:10){
  plot_tree(paste0("~/Documents/m_topo_",i,"/treemix_out"))
  plot_resid(paste0("~/Documents/m_topo_",i,"/treemix_out"), "~/Documents/Treemix_results/poporder_cat")
  
}
