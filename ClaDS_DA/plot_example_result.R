library(RPANDA)
library(coda)
load("ClaDS2_Cetacea.Rdata")


#the results are given as :
# the mcmc chains in coda_chain
# the phulogeny in tree
# the Maximum A Posteriori in MAPS[1:npar], where
#     MAPS[1] is the sigma parameter (stochasticity)
#     MAPS[2] is the alpha parameter (trend)
#     MAPS[3] is the epsilon parameter (turnover)
#     MAPS[4] is the lambda_0 parameter (initial speciation rate)
#     MAPS[5:npar] are the branch-specific speciation rates

plot_ClaDS_phylo(tree, MAPS[5:npar], log = T, lwd = 3)
title(main = paste0("sigma = ", round(MAPS[1],2), " ; alpha = ", round(MAPS[2],2),
                    " ; epsilon = ", round(MAPS[3],2), " ; lambda_0 = ", round(MAPS[4],2)))


#plot the chain for a given parameter id_par
id_par = 10
burn = 0
thin = 1
plot_chain = mcmc.list(lapply(1:3, function(i){
  mcmc(coda_chain[[i]][seq(max(1,floor(nrow(coda_chain[[1]])*burn)),nrow(coda_chain[[1]]), by = thin), id_par])
}))
plot(plot_chain)

#rmq : for rates (MAPS[5:npar]), the maps are computed on the log of the chains
plot_chain = mcmc.list(lapply(1:3, function(i){
  mcmc(log(coda_chain[[i]][seq(max(1,floor(nrow(coda_chain[[1]])*burn)),nrow(coda_chain[[1]]), by = thin), id_par]))
}))
plot(plot_chain)
