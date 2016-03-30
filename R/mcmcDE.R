library(parallel)

proposalGeneratorFactoryDE <- function(p=0.01,var=1e-6,burn=0,n.thin=0,decreasing.var=0,alpha=log(0.1)/1000,N.sigma=0){
  returnProposal <- function(chains,x,n){
    N=floor(n*(burn))
    npar=ncol(chains[[1]])-2
#     ind1.1=1
#     ind1.2=1
#     ind2.1=1
#     ind2.2=1
    
      
      ind1.1=sample(1:length(chains),1)
      ind2.1=sample(1:length(chains),1)
      if(n.thin==0 | (n-N)<=n.thin){
        # if(n.thin==0 | (n-N)<=n.thin){
        if(floor(N)==0){
#           ind1.2=n
#           ind2.2=n
          ind1.2=sample(1:n,1)
          ind2.2=sample(1:n,1)
        }else{
          ind1.2=sample(N:n,1)
          ind2.2=sample(N:n,1)
        }
      }else{
        ind1.2=sample(seq(max(1,floor(N)),n,length.out = n.thin),1)
        ind2.2=sample(seq(max(1,floor(N)),n,length.out = n.thin),1)
      }
      # }
    u=runif(1,0,1)
    
    if(u>p){
      gamma=2.36/sqrt(2*npar)
    }else{
      gamma=1
    }
    
    if(decreasing.var>0){
      e=rnorm(npar,0,var+decreasing.var*exp(alpha*n))
    }else{
      e=rnorm(npar,0,var)
      }
    
    
    rep=x+e+gamma*(chains[[ind1.1]][ind1.2,1:npar]-chains[[ind2.1]][ind2.2,1:npar])
    if(n<N.sigma){rep[1]=x[1]}
      
    return(rep)
  }
  
  return(list(returnProposal=returnProposal))
}

mcmcSamplerDE <- function(likelihood, Nchain=3, prior = NULL, startvalue, startmodel = NULL, iterations=10000, proposalGenerator = NULL, consoleupdates=1000000, thin = NULL){
  
  require(coda)
  require(compiler)
  require(MASS)
  
  ###############################################
  # Target definitions 
  
  if (is.null(prior)){
    prior <- function(x){
      return(0) 
    }
  }
  
  catchingLikelihood <- function(x){
    out <- tryCatch(
      {
        x = likelihood(x)
        if (x == Inf ){
          x = -Inf
          warning("Positive Inf values occured in the likelihood. Set to -Inf")
        }
        x 
      },
      error=function(cond){
        message("Problem in the likelihood")
        message(cond)
        return(-Inf)
      }
    )
    return(out)
  }
  
  catchingPrior <- function(x){
    out <- tryCatch(
      {
        prior(x)
      },
      error=function(cond) {
        message("Problem in the prior")
        message(cond)
        return(-Inf)
      }
    )
    if (out == Inf) out = -Inf
    return(out)
  }
  
  posterior <- function(x){
    priorResult = catchingPrior(x) # Checking first if outside the prior to save calculation time
    if (priorResult == -Inf) return(c(-Inf, -Inf))
    else{
      likelihoodResult = catchingLikelihood(x)
      return(c(likelihoodResult, likelihoodResult + priorResult))
    }
  }  
  
  posterior2 <- function(x){
    priorResult = catchingPrior(x) # Checking first if outside the prior to save calculation time
    if (priorResult == -Inf) return(-Inf)
    else return(catchingLikelihood(x) + priorResult)
  }
  
  numPars = length(startvalue[[1]]) + length(startmodel)
  
  if(is.null(proposalGenerator)){
    proposalGenerator = proposalGeneratorFactory(rep(1,numPars))
  }
  
  ####### CREATE CHAIN
  
  chains=list()
  currentLPs=list()
  for(i in 1:Nchain){
    chain = array(dim = c(1,numPars+2))
    chain[1,1:numPars] = c(startvalue[[i]], startmodel)
    colnames(chain) = c(1:numPars, "LL", "LP")
    chain[1, (numPars+1):(numPars+2)] = posterior(startvalue[[i]])
    currentLPs[[i]] = chain[1, (numPars+2)]
    chains[[i]]=chain
  }
  
  
  ##### Sampling
  
  classFields = list(
    likelihood = likelihood, 
    catchingLikelihood = catchingLikelihood,
    prior = prior,
    catchingPrior = catchingPrior,
    posterior = posterior,
    startvalue = startvalue, 
    numPars = numPars,
    indexLL = numPars + 1,
    indexLP = numPars + 2,
    currentLPs = currentLPs,
    chains = chains, 
    optimize = optimize, 
    consoleupdates=consoleupdates, 
    thin = thin,
    Nchain=Nchain,
    proposalGenerator = proposalGenerator,
    codaChain = NULL
  )
  
  class(classFields) <- append(class(classFields),"mcmcSampler")
  return(classFields)
}


getSamplesDE <- function(mcmcSampler, iterations, thin=NULL,nCPU=1){
  
  if(is.null(mcmcSampler$thin)){
    if(is.null(thin)){
      mcmcSampler$thin=1
    }else{
      mcmcSampler$thin=thin
    }
  }else{
    if(!is.null(thin)){
      if (!mcmcSampler$thin==thin){
        warning("thinning factor has been modified")
        mcmcSampler$thin=thin
      }
    }
  }
  lastvalue = nrow(mcmcSampler$chains[[1]])
  k=lastvalue
  currentchain=lapply(1:mcmcSampler$Nchain,function(i){mcmcSampler$chains[[i]][lastvalue,]})

  for (i in (lastvalue+1):(lastvalue+iterations/mcmcSampler$thin)){
    # for(j in 1:mcmcSampler$Nchain){    
    currentchain=mclapply(1:mcmcSampler$Nchain,function(j){
      current=currentchain[[j]]
      for (l in 1:mcmcSampler$thin){
      proposal = mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k)
      proposalEval <- mcmcSampler$posterior(proposal)
      probab = exp(proposalEval[2] - current[mcmcSampler$numPars+2])
      
      if (runif(1) < probab){
        current=c(proposal, proposalEval)
      }}
      return(current)
      },mc.cores = nCPU)
      
      # if( i %% mcmcSampler$thin == 0 ) {
        mcmcSampler$chains=lapply(1:mcmcSampler$Nchain,function(j){rbind(mcmcSampler$chains[[j]],currentchain[[j]])})
      flush.console()
      # }
      # }
    
    # if( i %% mcmcSampler$thin == 0 ) { 
      k=k+1
      # }
    if( (i-1)*mcmcSampler$thin %% mcmcSampler$consoleupdates == 0 ) cat("\r","MCMC in progress",(i-1)*mcmcSampler$thin,"of",iterations+mcmcSampler$thin*(lastvalue-1),"please wait!","\r")}
  
  message("Done.")
  mcmcSampler$codaChain = mcmc.list(lapply(1:mcmcSampler$Nchain,function(i){mcmc(mcmcSampler$chains[[i]])}))
  return(mcmcSampler)
}

burn.list=function(mcmcSampler,burn,thin=1,col=NULL,chain=NULL){
  rep=list()
  n=nrow(mcmcSampler$codaChain[[1]])
  if(is.null(chain)){chain=length(mcmcSampler$codaChain)}
  for (i in 1:chain){
    if(is.null(col)){rep[[i]]=mcmc(mcmcSampler$codaChain[[i]][seq(burn,n,thin),])
    }else{rep[[i]]=mcmc(mcmcSampler$codaChain[[i]][seq(burn,n,thin),col])}
  }
  return(mcmc.list(rep))
}

# ############test###################
# 
if (T){
  seed=92146169
  set.seed(seed)
  d=0
  a=birthdeath.tree.rateshift(1,0.1,d,sigma=0.03,taxa.stop = 100,condition="taxa",new_lamb_law = "normal",mu_min = d,mu_max=d,prune.extinct = T)
  tree=a$tree
  ntips=tree$Nnode+1
  nedges=2*tree$Nnode
  size_nt = 800
  
  true.rate=a$lamb$par[a$rates]
  f=1
  nCPU=1
  
  name_method="FFT" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nt,nt = size_nt,method=name_method , nCPU=nCPU) 
  
  start=list()
  N=1
  start=lapply(1:N,function(i){runif(n = nedges+3,0,0.01)})
  
  target <- function(x) {if(x[1]<1e-4){-Inf}else{fun(x[-(1:2)],x[1],x[2] ,f)}}
  target(start[[1]])
  testGenerator <- proposalGeneratorFactoryDE(var = 1e-6,p=0.1)
  sampler=mcmcSamplerDE(target,startvalue = start,proposalGenerator = testGenerator,Nchain = N,consoleupdates = 1)
  sampler=getSamplesDE(mcmcSampler = sampler,iterations = 200,thin = 1,nCPU=1)
  plot(sampler$codaChain)
}