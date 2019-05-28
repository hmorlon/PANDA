plot_ClaDS_chains = function(sampler,burn=1/2,thin=1,param=c("sigma","alpha","mu","LP")){
  nR=nrow(sampler$chains[[1]])
  npar=ncol(sampler$chains[[1]])-2
  rows =seq(max(1,floor(nR*(burn))),nR,by = thin)
  plot(mcmc.list(lapply(sampler$chains,function(x){mcmc(x[rows,param])})))
}