getMAPS_ClaDS0 = function(phylo, sampler, burn=1/2, thin=1){
  nR=nrow(sampler[[1]])
  npar=ncol(sampler[[1]])-3
  rows =seq(max(1,floor(nR*(burn))),nR,by = thin)
  
  chains=mcmc.list(lapply(sampler,function(x){mcmc(x[rows,])}))
  for(k in 1:length(chains)){
    for(l in 1:nrow(chains[[k]])){
      chains[[k]][l,4:npar]=get_rates(phylo,c(chains[[k]][l,3],
                                             chains[[k]][l,4:npar]+chains[[k]][l,2]))[-1]
    }
  }
  
  MAPS=sapply(1:npar, function(i){D=density(c(chains[[1]][,i],
                                              chains[[2]][,i],
                                              chains[[3]][,i]))
  return(D$x[which.max(D$y)])})
  MAPS[2:npar]=exp(MAPS[2:npar])
  return(MAPS)
}
