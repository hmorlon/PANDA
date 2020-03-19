library(RPANDA)

load("ClaDS2_Cetacea.Rdata")       # load the results in the environment


n_edges = nrow(tree$edge)                 # number of branches in the tree

sigma_inf = MAPS[1]                       # rate stochasticity
alpha_inf = MAPS[2]                       # trend in rate changes
epsilon_inf = MAPS[3]                     # turnover rate
l0_inf = MAPS[4]                          # initial speciation rate
rates_inf = MAPS[4+(1:n_edges)]           # branch-specific speciation rates

plot_ClaDS_phylo(tree,rates_inf)          # plot the tree
title(main = paste0("sigma", " = ", round(sigma_inf,2), " ;  ",
      "alpha", " = ", round(alpha_inf,2), " ;  ",
      "epsilon", " = ", round(epsilon_inf,2)))

plot_ClaDS0_chains(coda_chain, thin = 1, burn = 1/4, param = 1:4)
