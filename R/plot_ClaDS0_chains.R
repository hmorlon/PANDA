plot_ClaDS0_chains = function(sampler,burn=1/2,thin=1,param=c("sigma","alpha","l_0","LP")){
  plot(mcmc.list(lapply(sampler,function(x){mcmc(x[seq(max(1,floor(nrow(x)*(burn))),nrow(x),by = thin),param])})))
}