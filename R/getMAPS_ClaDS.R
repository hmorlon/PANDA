getMAPS_ClaDS = function(sampler, burn=1/2, thin=1){
  nR=nrow(sampler$chains[[1]])
  npar=ncol(sampler$chains[[1]])-2
  rows =seq(max(1,floor(nR*(burn))),nR,by = thin)
  
  chains=mcmc.list(lapply(sampler$chains,function(x){mcmc(x[rows,])}))
  
  for(k in 1:length(chains)){
    for(j in 1:nrow(chains[[1]])){
      chains[[k]][j,5:npar]=chains[[k]][j,4]+chains[[k]][j,2]*sampler$alpha_effect+
        chains[[k]][j,5:npar]
    }}
  
  MAPS=sapply(1:npar, function(i){D=density(c(chains[[1]][seq(1,nrow(chains[[1]]),10),i],
                                              chains[[2]][seq(1,nrow(chains[[1]]),10),i],
                                              chains[[3]][seq(1,nrow(chains[[1]]),10),i]))
  return(D$x[which.max(D$y)])})
  
  MAPS[-3]=exp(MAPS[-3])
  return(MAPS)
}
