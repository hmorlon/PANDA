scalar=function(x,y){
  sum(x*y)
}

norm2=function(x){
  sum(x^2)
}

proposalGeneratorFactoryDE_gibbs <- function(proba.gibbs,p=0.01,var=1e-6,burn=0,
                                             n.thin=0,decreasing.var=0,
                                             alpha=log(0.1)/200,N.sigma=0,allHyp=F){
  
  returnProposal <- function(chains,x,n,snook=F){
    N=n*(burn)
    npar=ncol(chains[[1]])-2
    if(length(var)==1){
      var=rep(var,npar)
    }
    if(length(decreasing.var)==1){
      decreasing.var=rep(decreasing.var,npar)
    }
    if(snook){
      ind1.1=sample(1:length(chains),1)
      ind2.1=sample(1:length(chains),1)
      ind3.1=sample(1:length(chains),1)
      if(n.thin==0 | (n-N)<=n.thin){
        
        if(floor(N)==0){
          
          ind1.2=sample(1:n,1)
          ind2.2=sample(1:n,1)
          ind3.2=sample(1:n,1)
        }else{
          ind1.2=sample(N:n,1)
          ind2.2=sample(N:n,1)
          ind3.2=sample(N:n,1)
        }
      }else{
        ind1.2=sample(seq(floor(N),n,length.out = n.thin),1)
        ind2.2=sample(seq(floor(N),n,length.out = n.thin),1)
        ind3.2=sample(seq(floor(N),n,length.out = n.thin),1)
      }
      # }
      
      
      if(any(decreasing.var>0)){
        e1=rnorm(npar,0,var+decreasing.var*exp(alpha*n))
        e2=0#rnorm(npar,0,var+decreasing.var*exp(alpha*n))
      }else{
        e1=rnorm(npar,0,var)
        e2=0#rnorm(npar,0,var)
      }
      
      z=chains[[ind3.1]][ind3.2,1:npar]-e1
      orth=x-z
      if(all(orth==0)){
        snook=orth
      }else{
        orth=orth/sqrt(norm2(orth))
        snook=(scalar(orth,chains[[ind1.1]][ind1.2,1:npar]-chains[[ind2.1]][ind2.2,1:npar])-e2)*orth
      }
      
      
      gamma=(2.36/sqrt(2))
      
      
      rep=x
      # print(x)
      rep=rep+(gamma*snook)
      if(n<N.sigma){rep[1]=x[1]}
      
      return(list(prop=rep,z=z))
      
    }else{
        gibbs=sample(1:npar,sample(1:npar,1,prob = proba.gibbs))
        if(allHyp){
          if(1 %in% gibbs | 2 %in% gibbs | 3 %in% gibbs){
          gibbs=1:npar
        }}
        
        ind1.1=sample(1:length(chains),1)
        ind2.1=sample(1:length(chains),1)
        if(n.thin==0 | (n-N)<=n.thin){
          
          if(floor(N)==0){
            
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
        # gamma=2.36/sqrt(2*npar)
        # gamma=runif(1,0,2*(2.36/sqrt(2*length(gibbs))))
        gamma=rexp(1,1/(2.36/sqrt(2*length(gibbs))))
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
  }
  
  return(list(returnProposal=returnProposal))
}


mcmcSamplerDE_gibbs <- function(likelihood,proba.gibbs, Nchain=3, 
                                startvalue, startmodel = NULL, 
                                iterations=10000, proposalGenerator = NULL, consoleupdates=1000000, thin = NULL){
  

  
  ###############################################
  # Target definitions 

   prior <- function(x){
      return(0) 
  }

  
  
  numPars = length(startvalue[[1]]) + length(startmodel)
  
  if(is.null(proposalGenerator)){
    proposalGenerator = proposalGeneratorFactoryDE_gibbs(rep(1,numPars))
  }
  
  ####### CREATE CHAIN
  
  chains=list()
  currentLPs=list()
  former=list()
  for(i in 1:Nchain){
    if(i==1 | !identical(startvalue[[i]],startvalue[[1]])){
      chain = array(dim = c(1,numPars+2))
      chain[1,1:numPars] = c(startvalue[[i]], startmodel)
      colnames(chain) = c(1:numPars, "LL", "LP")
      form=likelihood(startvalue[[i]])
      chain[1, (numPars+1):(numPars+2)] = c(form$LP,form$LP)
    }else{
      chain=chains[[1]]
      form=former[[1]]
    }
    
    former[[i]]=form
    currentLPs[[i]] = chain[1, (numPars+2)]
    chains[[i]]=chain
  }
  
  
  ##### Sampling
  
  classFields = list(
    likelihood = likelihood, 
    startvalue = startvalue, 
    numPars = numPars,
    currentLPs = currentLPs,
    chains = chains, 
    consoleupdates=consoleupdates, 
    thin = thin,
    Nchain=Nchain,
    proposalGenerator = proposalGenerator,
    post=likelihood,
    former=former
  )
  
  class(classFields) <- append(class(classFields),"mcmcSampler")
  return(classFields)
}

prepare_ClaDS=function(tree,sample_fraction, Nchain=3,model_id="ClaDS2", res_ClaDS0=NULL, l0 = 0.1, s0 = 1, nlambda = 1000, nt = 30){
  nedges=nrow(tree$edge)
  min_sig_Cl0=0
  # the likelihood and posterior functions
  if(model_id == "ClaDS1"){
    likelihood = createLikelihood_ClaDS1(tree,nlambda = nlambda,nt = nt, conv=1e-5) # the likelihood function
  }else{ 
    likelihood = createLikelihood_ClaDS2(tree,nlambda = nlambda,nt = nt, conv=1e-5) # the likelihood function
  }
  relToAbs = likelihood$relToAbs
  likelihood = likelihood$ll
  prior=function(x){0}  # the prior function, here a flat prior
  
  alpha_effect=relToAbs(c(0,rep(1,tree$Nnode*2)))[-1]
  
  posterior=function(param,former=NULL, ae = alpha_effect){        # the posterior function
    param2=param
    param2[1]=exp(param2[1])
    param2[2]=exp(param2[2])
    param2[-(1:4)]=param2[4]+param[-(1:4)]+ae*param[2]
    if(param2[1]<2){t=try(likelihood(param2,sample_fraction,former),silent = T)
    if(inherits(t,"try-error")){
      return(list(LP=-Inf))
      
    }else if(is.nan(t$LL) | is.infinite(t$LL)){
      return(list(LP=-Inf))
    }else{
      Pr=prior(param2)
      t$Pr=Pr
      t$LP=t$LL+Pr
      return(t)}}else{
        return(list(LP=-Inf))
        
      }}
  
  if(is.null(res_ClaDS0)){
    start=lapply(1:Nchain,function(i){c(sigma=log(s0),alpha=log(1),mu=0,log(l0),lambda=c(rnorm(nedges,0,0))) })#rel.to.abs(tree,c(log(0.001),rep(log(A[i]),nedges)))) })
    decreasing.var=c(rep(1e-3,2),rep(1e-1,2),rep(1e-2,nedges))
  }else{
      start=lapply(1:Nchain,function(i){x=res_ClaDS0[[i]]$chain[nrow(res_ClaDS0[[i]]$chain),1:(nedges+3)]
      return(c(log(max(min_sig_Cl0,x[1])),x[2],runif(1),relToAbs(x[-(1:2)])))})
      decreasing.var=sapply(1:3,function(i){res_ClaDS0[[i]]$finetune})
      decreasing.var=rowMeans(decreasing.var)
      decreasing.var=0.1*c(decreasing.var[1:2],1,decreasing.var[-(1:2)])
  }
  
  npar=length(start[[1]])
  g=3 #mean number of parameter updated at each iteration
  G=exp(-(1:(npar))/g)
  testGenerator <- proposalGeneratorFactoryDE_gibbs(decreasing.var=decreasing.var,alpha=log(0.1)/200,
                                                    var=c(1e-4,1e-4,rep(1e-4,2),rep(1e-4,nedges)),
                                                    proba.gibbs=c(G))
  
  sampler=mcmcSamplerDE_gibbs(posterior,startvalue = start,proposalGenerator = testGenerator,Nchain = Nchain,consoleupdates = 1)
  sampler$alpha_effect=alpha_effect
  sampler$relToAbs=relToAbs
  return(sampler)
}

add_iter_ClaDS <- function(mcmcSampler, iterations, thin=NULL,nCPU=1){
  snookProb=0.05
  colnames=c("sigma","alpha","mu","l_0")
  post=mcmcSampler$post
  alpha_effect=mcmcSampler$alpha_effect
  relToAbs=mcmcSampler$relToAbs
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
    currentchain=mclapply(1:mcmcSampler$Nchain,function(j){
      current=mcmcSampler$chains[[j]][k,]
      form=former[[j]]
      for (l in 1:mcmcSampler$thin){
        u=runif(1,0,1)
        if(u < snookProb){
          proposal = try(mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k,snook=T))
          proposalEval <- try(post(proposal$prop,form, ae = mcmcSampler$alpha_effect),silent = T)
          if(!is.null(proposal) & !inherits(proposal,"try-error") & !inherits(proposalEval,"try-error") & inherits(proposalEval,"list")){
            if(!is.nan(proposalEval$LP) & proposalEval$LP<Inf & proposalEval$LP>-Inf){
              probab = exp(proposalEval$LP + ((mcmcSampler$numPars-1)/2)*(log(norm2(proposal$prop-proposal$z))-
                                                                            log(norm2(current[1:mcmcSampler$numPars]-proposal$z))) - 
                             current[mcmcSampler$numPars+2])
              
              if ( runif(1) < probab){
                current=c(proposal$prop, proposalEval$LP,proposalEval$LP)
                form=proposalEval
              }}}
        }else{
          proposal = try(mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k,snook=F))
          proposalEval <- try(post(proposal,form),silent = T)
          if(!is.null(proposal) & !inherits(proposal,"try-error") & !inherits(proposalEval,"try-error") & inherits(proposalEval,"list")){
            if(!is.nan(proposalEval$LP) & proposalEval$LP<Inf & proposalEval$LP>-Inf){
              probab = exp(proposalEval$LP - current[mcmcSampler$numPars+2])
              
              if ( runif(1) < probab){
                current=c(proposal, proposalEval$LP,proposalEval$LP)
                form=proposalEval
              }}}
        }
      }
      return(list(current=current, former=form))
    },mc.cores = nCPU)
    
    if(sum(!sapply(currentchain,is.null))==mcmcSampler$Nchain ) {
      mcmcSampler$chains=lapply(1:mcmcSampler$Nchain,function(j){rbind(mcmcSampler$chains[[j]],currentchain[[j]]$current)})
      mcmcSampler$former=lapply(1:mcmcSampler$Nchain,function(j){currentchain[[j]]$former})
      
      flush.console()
      k=k+1
      i=i+1
    }
    if( i %% mcmcSampler$consoleupdates == 0 ) cat("\r","MCMC in progress",(i-2)*mcmcSampler$thin,"of",iterations+mcmcSampler$thin*(lastvalue-1),"please wait!","\r")}
  
  message("Done.")
  for(i in 1:mcmcSampler$Nchain){
    colnames(mcmcSampler$chains[[i]])[1:length(colnames)]=colnames
  }
  mcmcSampler$chains = mcmc.list(lapply(1:mcmcSampler$Nchain,function(i){mcmc(mcmcSampler$chains[[i]])}))
  mcmcSampler$alpha_effect=alpha_effect
  mcmcSampler$relToAbs=relToAbs
  return(mcmcSampler)
}

run_ClaDS = function(tree,sample_fraction,iterations, thin = 50, file_name = NULL, it_save = 1000,
                     model_id="ClaDS2", nCPU = 1, mcmcSampler = NULL, ...){
  args = list(...)
  
  if (is.null(args$Nchain)) args$Nchain=3
  if (is.null(args$l0)) args$l0=0.1
  if (is.null(args$s0)) args$s0=1
  if (is.null(args$nlambda)) args$nlambda=1000
  if (is.null(args$nt)) args$nt=30
  
  if (is.null(mcmcSampler)){
    mcmcSampler = prepare_ClaDS(tree=tree,          
                                sample_fraction=sample_fraction,
                                Nchain = args$Nchain,        
                                nlambda = args$nlambda,    
                                nt = args$nt,            
                                model_id=model_id,  
                                res_ClaDS0 = args$sampler_ClaDS0)
  }

  n_run = floor(iterations / it_save)
  
  if(n_run > 0){
    for (i in 1:n_run){
      mcmcSampler = add_iter_ClaDS(           
        mcmcSampler = mcmcSampler,     
        iterations = it_save,            
        thin = thin,  
        nCPU = nCPU)   
      
      if (! is.null(file_name)){
        save(mcmcSampler, file = file_name)
      }
      
      #print(paste0("iteration ", i * it_save," out of ", iterations))
    }
  }
  
  if((n_run * it_save)< iterations){
    mcmcSampler = add_iter_ClaDS(           
      mcmcSampler = mcmcSampler,     
      iterations = iterations - n_run *it_save,            
      thin = thin,  
      nCPU = nCPU)   
    
    if (! is.null(file_name)){
      save(mcmcSampler, file = file_name)
    }
    
    #print(paste0("iteration ", iterations," out of ", iterations))
  }

  return(mcmcSampler)
  
}
