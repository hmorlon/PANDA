library(parallel)

proposalGeneratorFactoryDE_gibbs <- function(proba.gibbs,p=0.01,var=1e-6,burn=0,n.thin=0,decreasing.var=0,alpha=log(0.1)/1000,N.sigma=0){

  returnProposal <- function(chains,x,n){
    N=n*(burn)
    npar=ncol(chains[[1]])-2
    if(length(var)==1){
      var=rep(var,npar)
    }
    if(length(decreasing.var)==1){
      decreasing.var=rep(decreasing.var,npar)
    }
    gibbs=sample(1:npar,sample(1:npar,1,prob = proba.gibbs))
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
        ind1.2=sample(seq(floor(N),n,length.out = n.thin),1)
        ind2.2=sample(seq(floor(N),n,length.out = n.thin),1)
      }
      # }
    u=runif(1,0,1)
    
    if(u>p){
      gamma=2.36/sqrt(2*length(gibbs))
    }else{
      gamma=1
    }
    
    if(any(decreasing.var[gibbs]>0)){
      e=rnorm(length(gibbs),0,var[gibbs]+decreasing.var[gibbs]*exp(alpha*n))
    }else{
      e=rnorm(length(gibbs),0,var[gibbs])
      }
    
    
    rep=x
    # print(x)
    rep[gibbs]=rep[gibbs]+e+(gamma*(chains[[ind1.1]][ind1.2,1:npar]-chains[[ind2.1]][ind2.2,1:npar]))[gibbs]
    if(n<N.sigma){rep[1]=x[1]}
      
    return(rep)
  }
  
  return(list(returnProposal=returnProposal))
}

mcmcSamplerDE_gibbs <- function(likelihood,proba.gibbs, Nchain=3, prior = NULL, startvalue, startmodel = NULL, iterations=10000, proposalGenerator = NULL, consoleupdates=1000000, thin = NULL){
  
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
  
#   catchingLikelihood <- function(x){
#     out <- tryCatch(
#       {
#         x = likelihood(x)
#         if (x == Inf ){
#           x = -Inf
#           warning("Positive Inf values occured in the likelihood. Set to -Inf")
#         }
#         x 
#       },
#       error=function(cond){
#         message("Problem in the likelihood")
#         message(cond)
#         return(-Inf)
#       }
#     )
#     return(out)
#   }
#   
#   catchingPrior <- function(x){
#     out <- tryCatch(
#       {
#         prior(x)
#       },
#       error=function(cond) {
#         message("Problem in the prior")
#         message(cond)
#         return(-Inf)
#       }
#     )
#     if (out == Inf) out = -Inf
#     return(out)
#   }
#   
#   posterior <- function(x){
#     priorResult = catchingPrior(x) # Checking first if outside the prior to save calculation time
#     if (priorResult == -Inf) return(c(-Inf, -Inf))
#     else{
#       likelihoodResult = catchingLikelihood(x)
#       return(c(likelihoodResult, likelihoodResult + priorResult))
#     }
#   }  
#   
#   posterior2 <- function(x){
#     priorResult = catchingPrior(x) # Checking first if outside the prior to save calculation time
#     if (priorResult == -Inf) return(-Inf)
#     else return(catchingLikelihood(x) + priorResult)
#   }
  
  numPars = length(startvalue[[1]]) + length(startmodel)
  
  if(is.null(proposalGenerator)){
    proposalGenerator = proposalGeneratorFactory(rep(1,numPars))
  }
  
  ####### CREATE CHAIN
  
  chains=list()
  currentLPs=list()
  former=list()
  for(i in 1:Nchain){
    chain = array(dim = c(1,numPars+2))
    chain[1,1:numPars] = c(startvalue[[i]], startmodel)
    colnames(chain) = c(1:numPars, "LL", "LP")
    form=likelihood(startvalue[[i]])
    chain[1, (numPars+1):(numPars+2)] = c(form$logLik,form$logLik)
    former[[i]]=form
    currentLPs[[i]] = chain[1, (numPars+2)]
    chains[[i]]=chain
  }
  
  
  ##### Sampling
  
  classFields = list(
    likelihood = likelihood, 
    # catchingLikelihood = catchingLikelihood,
    prior = prior,
#     catchingPrior = catchingPrior,
#     posterior = posterior,
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
    codaChain = NULL,
    post=likelihood,
    former=former
  )
  
  class(classFields) <- append(class(classFields),"mcmcSampler")
  return(classFields)
}

getSamplesDE_gibbs <- function(mcmcSampler, iterations, post, thin=NULL,nCPU=1,colnames=c("sigma","mu","l_0")){
  post=mcmcSampler$post
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
  former=mcmcSampler$former
  i=lastvalue+1

  while (i <((lastvalue+iterations/mcmcSampler$thin)+1)){
    # for(j in 1:mcmcSampler$Nchain){    
    currentchain=mclapply(1:mcmcSampler$Nchain,function(j){
      current=mcmcSampler$chains[[j]][k,]
      form=former[[j]]
      for (l in 1:mcmcSampler$thin){
      proposal = try(mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k))
      proposalEval <- try(post(proposal,form),silent = T)
      if(!is.null(proposal) & !inherits(proposal,"try-error") & !inherits(proposalEval,"try-error") & inherits(proposalEval,"list")){
#         if(is.nan(proposalEval$logLik)){debug(Phi)
#           post(proposal,form)
#           print(proposal)}
        if(!is.nan(proposalEval$logLik) & proposalEval$logLik<Inf & proposalEval$logLik>-Inf){
      # print(proposalEval$logLik)
      probab = exp(proposalEval$logLik - current[mcmcSampler$numPars+2])
      
      if (is.nan(probab) | runif(1) < probab){
        current=c(proposal, proposalEval$logLik,proposalEval$logLik)
        form=proposalEval
      }}}
      }
      return(list(current=current, former=form))
      },mc.cores = nCPU)
      
      if(sum(!sapply(currentchain,is.null))==mcmcSampler$Nchain ) {
        mcmcSampler$chains=lapply(1:mcmcSampler$Nchain,function(j){rbind(mcmcSampler$chains[[j]],currentchain[[j]]$current)})
        mcmcSampler$former=lapply(1:mcmcSampler$Nchain,function(j){currentchain[[j]]$former})
        # print(currentchain[[1]]$former$logLik)
      flush.console()
      # }
      # }
    
    # if( i %% mcmcSampler$thin == 0 ) { 
      k=k+1
      i=i+1
      }
    if( (i-1)*mcmcSampler$thin %% mcmcSampler$consoleupdates == 0 ) cat("\r","MCMC in progress",(i-1)*mcmcSampler$thin,"of",iterations+mcmcSampler$thin*(lastvalue-1),"please wait!","\r")}
  
  message("Done.")
  for(i in 1:mcmcSampler$Nchain){
    colnames(mcmcSampler$chains[[i]])[1:length(colnames)]=colnames
  }
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

createLikelihood_death_gibbs <- function(phylo, root_depth=0, relative = F,nt=1000,M=2,method="expoRkit",banded=T,nlambda_max=5000){
  
  nbtips = Ntip(phylo)
  edge = phylo$edge
  edge.length = phylo$edge.length
  roots=c()
  
  type = rep(NA, Ntip(phylo))
  parents = rep(NA, Ntip(phylo))
  offspring=matrix(rep(NA, 2*(2*Ntip(phylo)-2)),nrow=2)
  
  for (i in 1:(2*(nbtips - 1))) {
    parent=which(phylo$edge[,2]==phylo$edge[i,1])
    if(length(parent)==0){
      roots=c(roots,i)
      parents[i] = 0
    }else{
      parents[i] = parent
      if(is.na(offspring[1,parent])){offspring[1,parent]=i}else{offspring[2,parent]=i}
    } 
    if (!(phylo$edge[i,2]<=nbtips)){
      type[i] = 2
    }else{
      type[i] = 1       
    }
  }
  
  ancestors=list()
  for(i in 1:nedges){
    if(i %in% roots){
      ancestors[[i]]=-1
    }else{
      ancestors[[i]]=c(parents[i],ancestors[[parents[i]]][ancestors[[parents[i]]]>0])
    }
  }
  
  nodeprof=node.depth.edgelength(phylo)
  nodeprof=max(nodeprof)-nodeprof+root_depth
  tf=max(nodeprof)
  
  relToAbs <- function(lambda){
    lambda2=lambda
    for(i in 2:length(lambda)){
      lambda2[i]=lambda2[parents[i-1]+1]+lambda2[i]
    }
    return(lambda2)
  }
  nodeprof=node.depth.edgelength(phylo)
  nodeprof=max(nodeprof)-nodeprof+root_depth
  tf=max(nodeprof)
  
  
  ll<-function(lambda, sigma, mu,f,former=NULL,nlambda=min(nlambda_max,5*M/(sigma))){
    
    lambda2=lambda
    lambda2=relToAbs(lambda2)
    
    if (any(lambda2 <= 0) | any(lambda2>M)) return(list(lambda=lambda2,mu=mu,sigma=sigma,phi=NULL,PsiKhi=NULL,f=f,logLik=-Inf))
    if (any(c(mu,f,sigma) < 0)) return(list(lambda=lambda2,mu=mu,sigma=sigma,phi=NULL,PsiKhi=NULL,f=f,logLik=-Inf))
    eval=is.null(former)
    if(!eval){eval=((former$sigma != sigma)|(former$mu != mu)|(former$f != f))}
    
    if(eval){
      PsiKhi=rep(0,nedges+1)
      if(method=="FFT"){
        phi=Phi_FFT(sigma,M,nlambda,mu,f,tf=tf,by=tf/nt)
      }else{
        phi=Phi(sigma,M,nlambda,mu,f,tf=tf,by=tf/nt)
      }
      for(i in 1:(length(lambda2)-1)){
        if(type[i]==1){
          m=Khi(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,banded = banded)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          PsiKhi[i+1]=log(m[ind])
        }else{
          m=Khi(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                lambda2 = lambda2[offspring[2,i]+1],func = "Zeta",nt=nt,method=method,banded = banded)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          PsiKhi[i+1]=log(m[ind])}}
      if(root_depth==0){
        PsiKhi[1]=sum(dnorm(lambda2[roots+1],mean=lambda2[1],sd=sigma,log=T))
      }else{
        m=Khi(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
              lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt)
        ind=which.min(abs(lambda2[1]-phi$lambda))
        PsiKhi[1]=log(m[ind])
      }
    }else{
      change=which((lambda2-former$lambda)!=0)-1
      change=unique(c(change,sapply(change, function(x){if(x>0){parents[x]}else{-1}})))
      change=change[change>(-1)]
      PsiKhi=former$PsiKhi
      phi=former$phi
      for( i in change){
        if(i ==0){
          if(root_depth==0){
            PsiKhi[1]=sum(dnorm(lambda2[roots+1],mean=lambda2[1],sd=sigma,log=T))
          }else{
            m=Khi(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
                  lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,method=method,banded = banded)
            ind=which.min(abs(lambda2[1]-phi$lambda))
            PsiKhi[1]=log(m[ind])
          }
        }else{
          if(type[i]==1){
            m=Khi(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,banded = banded)
            ind=which.min(abs(lambda2[i+1]-phi$lambda))
            PsiKhi[i+1]=log(m[ind])
          }else{
            m=Khi(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                  lambda2 = lambda2[offspring[2,i]+1],func = "Zeta",nt=nt,method=method,banded = banded)
            ind=which.min(abs(lambda2[i+1]-phi$lambda))
            PsiKhi[i+1]=log(m[ind])
          }
        }
      }
    }
    rep=list(lambda=lambda2,mu=mu,sigma=sigma,phi=phi,PsiKhi=PsiKhi,f=f,logLik=sum(PsiKhi))
    return(rep)
  }
  
  return(ll)
  
}

# ############test###################

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
  fun=createLikelihood_death2(tree,relative = T,nlambda_max = 1000,nt = 1000,M=2,banded = 10000,method = "FFT")
  target <- function(x,former=NULL) {if(x[1]<1e-4){-Inf}else{fun(x[-(1:2)],x[1],x[2] ,1,former=former)}}
  
  start=list()
  N=1
  start=lapply(1:N,function(i){runif(n = nedges+3,0,0.01)})
  
  target(start[[1]])
  testGenerator <- proposalGeneratorFactoryDE_gibbs(var = 1e-6,p=0.1)
  sampler=mcmcSamplerDE_gibbs(target,startvalue = start,proposalGenerator = testGenerator,Nchain = N,consoleupdates = 1)
  sampler=getSamplesDE_gibbs(mcmcSampler = sampler,iterations = 200,thin = 1,nCPU=1)
  
  g=10
  testGenerator <- proposalGeneratorFactoryDE_gibbs(var=1e-5,proba.gibbs=c(rep(1,g),rep(0,nedges+3-g)),p=0.1)
  sampler=mcmcSamplerDE_gibbs(target,startvalue = start,proposalGenerator = testGenerator,Nchain = N,consoleupdates = 1)
  sampler=getSamplesDE_gibbs(mcmcSampler = sampler,iterations = 200,thin = 50, nCPU = N)
  plot(sampler$codaChain)
}