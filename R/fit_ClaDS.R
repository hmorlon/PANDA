#### general internal functions ####


scalar=function(x,y){
  sum(x*y)
}

norm2=function(x){
  sum(x^2)
}

proposalGeneratorFactoryDE_gibbs <- function(proba.gibbs,p=0.01,var=1e-6,burn=0,
                                             n.thin=0,decreasing.var=0,
                                             alpha=log(0.1)/200,N.sigma=0,allHyp=FALSE){
  
  returnProposal <- function(chains,x,n,snook=FALSE){
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
    if(param2[1]<2){t=try(likelihood(param2,sample_fraction,former),silent = TRUE)
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

add_iter_ClaDS <- function(mcmcSampler, iterations, thin=NULL,nCPU=1, verbose=TRUE){
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
          proposal = try(mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k,snook=TRUE))
          proposalEval <- try(post(proposal$prop,form, ae = mcmcSampler$alpha_effect),silent = TRUE)
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
          proposal = try(mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k,snook=FALSE))
          proposalEval <- try(post(proposal,form),silent = TRUE)
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
    if( i %% mcmcSampler$consoleupdates == 0 ) if (verbose) cat("\r","MCMC in progress",(i-2)*mcmcSampler$thin,"of",iterations+mcmcSampler$thin*(lastvalue-1),"please wait!","\r")}
  
  message("Done.")
  for(i in 1:mcmcSampler$Nchain){
    colnames(mcmcSampler$chains[[i]])[1:length(colnames)]=colnames
  }
  mcmcSampler$chains = mcmc.list(lapply(1:mcmcSampler$Nchain,function(i){mcmc(mcmcSampler$chains[[i]])}))
  mcmcSampler$alpha_effect=alpha_effect
  mcmcSampler$relToAbs=relToAbs
  return(mcmcSampler)
}


InterpolatedPhi=function(phi,t){
  tind=which.min(abs(phi$fun[,1]-t))
  if (t >= phi$fun[nrow(phi$fun),1]) {
    return(phi$fun[nrow(phi$fun),-1])
  } else if (t <= phi$fun[1,1]) {
    return(phi$fun[1,-1])
  } else if (t >= phi$fun[tind,1]) {
    alpha = (t - phi$fun[tind,1]) / (phi$fun[tind+1,1] - phi$fun[tind,1])
    return(alpha*phi$fun[tind+1,-1] + (1-alpha)*phi$fun[tind,-1])
  } else {
    alpha = (t - phi$fun[tind-1,1]) / (phi$fun[tind,1] - phi$fun[tind-1,1])
    return(alpha*phi$fun[tind,-1] + (1-alpha)*phi$fun[tind-1,-1])
  }
}


InterpolatedKhi=function(phi,m,lambda){
  ind=which.min(abs(lambda-phi$lambda))
  
  if (lambda >= phi$lambda[ind]) {
    if(ind==(length(phi$lambda))){
      return(m[ind])
    }
    alpha = (lambda - phi$lambda[ind]) / (phi$lambda[ind+1] - phi$lambda[ind])
    return(log(alpha*(m[ind+1]) + (1-alpha)*(m[ind])))
  } else {
    if(ind==1){
      return(m[ind])
    }
    alpha = (lambda - phi$lambda[ind-1]) / (phi$lambda[ind] - phi$lambda[ind-1])
    return(log(alpha*(m[ind]) + (1-alpha)*(m[ind-1])))
  }
}

abs_to_rel=function(tree,true.rate,lambda_0){
  rate=c(lambda_0,rep(0,length(true.rate)))
  ntips=tree$Nnode+1
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,1]==ntips+1){
      rate[i+1]=true.rate[i]-lambda_0
    }else{
      parent=which(tree$edge[,2]==tree$edge[i,1])
      rate[i+1]=true.rate[i]-true.rate[parent]
    }
  }
  return(rate)
}


#### for ClaDS1 ###

Phi_ClaDS1=function(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1){
  lambdaIs = seq(log(mlambda),log(Mlambda),length=nlambda+1)
  
  firstRow = -(lambdaIs-log(mlambda))+log(alpha)
  firstCol = lambdaIs-log(mlambda)+log(alpha)
  if(sigma==0){
    row_indices = which.min(abs(firstRow))
    firstRow = rep(0,nlambda+1)
    firstRow[row_indices] = 1
    col_indices = which.min(abs(firstCol))
    firstCol = rep(0,nlambda+1)
    firstCol[col_indices] = 1
  }else{
    firstRow = exp(-firstRow^2/(2*sigma^2))
    firstCol = exp(-firstCol^2/(2*sigma^2))
  }
  M = MPhiFFT(firstCol, firstRow, method="FFT")
  
  ini=rep((1-f),nlambda+1)
  #   ini[1]=0
  #   ini[length(ini)]=0
  ind=(1:length(lambdaIs))
  expLambda=exp(lambdaIs)
  
  MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
  dPhi=function(t,y,parms){
    
    dy=expLambda*(MATVECT(M, y)^2-y)+mu*(1-y)
    return(list(dy))
  }
  
  times <- seq(from = tini, to =tf, by = by)
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL, method="daspk", atol=1e-3)
  out[,-1][out[,-1]>(1-1e-3)]=1
  out[,-1][out[,-1]<0]=0
  
  return(list(lambda=lambdaIs,expLambda=expLambda,fun=out,mu=mu,f=f,sigma=sigma,M=M,func="Phi",alpha=alpha,ind=ind))
}

MagnusExpansion_ClaDS1=function(phi,ini,tini,tend,step,conv=1e-10,mask_val=-20,timeIni,timeEnd){
  expLambda=phi$expLambda
  mu=phi$mu
  timePhi=phi$fun[,1]
  M=phi$M
  iniRep=ini
  if(tend<(tini+step)) tend=tini+step
  MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
  tStep=step*diff(timePhi)[1]
  S=seq(timeIni, max(timeIni,timeEnd-tStep), by=tStep)
  if(length(S)==1) {
    h=timeEnd-timeIni
  }else{
    h = S[2]-S[1]
  }
  D0 = expLambda+mu
  mask = - (expLambda+mu)*h < mask_val
  D3 = expLambda*MATVECT(M, phi$fun[tini+step,-1])
  for (i in S) {
    # Compute iteratively exp(Omega)*V
    
    # Magnus expansion of order 4 with equispaced points:
    #    B1 = B(Tn), B2 = B(Tn+h/2), B3 = B(Tn+h)
    #    Omega(h) = (h/6)*(B1 + 4*B2 + B3) - (h*h/12) [B1,B3]
    #    Y(n+1) = exp(Omega(h)) Y(n)
    #
    # We compute Y(n+1) = exp(Omega(h)) %*% Y(n) by an iterative process
    #    OV(k+1) = Omega %*% OV(k) / (k+1)
    #    EXPOV(k+1) = EXPOV(k) + OV(k+1)
    # Series EXPOV(k) converges towards exp(Omega(h)) %*% Y(0), and we only need
    # to explicit the product of Omega(h) by a vector.
    #
    # Here, B(t) = -diag(expLambda+mu) + (2*expLambda*diag(M %*% Phi(t))) * M
    #    D0 = diag(expLambda+mu)
    #    D1 = expLambda*diag(M %*% Phi(Tn)),
    #    D2 = expLambda*diag(M %*% Phi(Tn+h/2))
    #    D3 = expLambda*diag(M %*% Phi(Tn+h))
    #    D123 = D1+4*D2+D3 = expLambda*diag(M %*% (Phi(Tn) + 4 Phi(Tn+h/2) + Phi(Tn+h)))
    #
    # This gives
    #    Omega(h)  = - h*D0 + (h/3)*D123 * M - (h*h/12) [B1,B3]
    #
    # First two terms are straightforward, we now have to explicit the last one  [B1,B3] = B1*B3 - B3*B1
    # Recall that
    #    B1 = - D0 + 2*D1*M
    #    B3 = - D0 + 2*D3*M
    # Thus,
    #    B1 * B3 = D0^2 - 2*D0*D3*M - 2*D1*M*D0 + 4*D1*M*D3*M
    #    B3 * B1 = D0^2 - 2*D0*D1*M - 2*D3*M*D0 + 4*D3*M*D1*M
    #    [B1,B3] = - 2*D0*(D3-D1)*M + 2*(D3-D1)*M*D0 + 4*(D1*M*D3*M - D3*M*D1*M)
    #            = - 2*(D3-D1)*(D0*M - M*D0) + 4*(D1*M*D3*M - D3*M*D1*M)
    # Eventually,
    #    Omega(h) = - h*D0 + (h/3)*D123*M - (h*h/12)*[B1,B3]
    #             = - h*D0 + (h/3)*D123*M + (h*h/6)*(D3-D1)*(D0*M - M*D0) - (h*h/3)*(D1*M*D3*M - D3*M*D1*M)
    
    if(length(S)==1){
      phi123 = InterpolatedPhi(phi,timeIni) + 4 * InterpolatedPhi(phi,(timeEnd+timeIni)/2) + InterpolatedPhi(phi,timeEnd)
      D1 = D3
      D3 = expLambda*MATVECT(M, InterpolatedPhi(phi,timeEnd))
    }else{
      phi123 = InterpolatedPhi(phi,i)+ 4 * InterpolatedPhi(phi,i+tStep/2) + InterpolatedPhi(phi,i+tStep)
      D1 = D3
      D3 = expLambda*MATVECT(M,InterpolatedPhi(phi,i+tStep))
    }
    D123 = expLambda*MATVECT(M,phi123)
    
    ini[mask] = 0
    EXPOV=ini
    OV=ini
    k=1
    epsnormv = conv * sum(abs(OV))
    end=0
    while (sum(abs(OV)) > epsnormv & end < 1e3){
      end=end+1
      MOV = MATVECT(M,OV)
      OV = - h*D0*OV + (h/3)*D123*MOV + (h*h/6)*(D3-D1)*(D0*MOV - MATVECT(M,D0*OV)) - (h*h/3)*(D1*MATVECT(M,D3*MOV)-D3*MATVECT(M,D1*MOV))
      OV = OV/k
      OV[mask] = 0
      EXPOV = EXPOV + OV
      k=k+1
    }
    if(end==1e3){EXPOV[]=0}
    EXPOV[mask] = 0
    ini = EXPOV
  }
  EXPOV[mask] = exp(- (expLambda[mask]+mu)*(timeEnd - timeIni))*iniRep[mask]
  return(EXPOV)
}

Khi_ClaDS1=function(phi,s,t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
                    timePhi=phi$fun[,1],nt=1000,method="Higham08.b",banded=0,sparse=FALSE,interval=2,conv=1e-10){
  
  if((0.99*t)>phi$fun[nrow(phi$fun),1]) return(rep(1e-300,length(phi$expLambda)))
  if( !(func %in% c("Khi","Psi","Ksi","Zeta"))){stop("the function is not known")}
  if(any(is.na(phi$fun))){stop("Phi is NA")}
  
  tini=which.min(abs(timePhi-s))
  tend=which.min(abs(timePhi-t))
  nlambda=length(lambdas)
  expLambda=phi$expLambda
  
  if(func=="Khi" | func=="Psi"){
    
    ini=1-phi$fun[tini,-1]  
    
  }else{
    #on transforme la proba d'être dans un intervalle (de taille lambdas[2]-lambdas[1]) en densité 
    ini=expLambda*dnorm(lambda1,lambdas+log(phi$alpha),phi$sigma)*dnorm(lambda2,lambdas+log(phi$alpha),phi$sigma)
    
  }
  
  if(method == "ode"){
    MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
    dKhi=function(t,y,parms){
      tindex=which.min(abs(timePhi-t))
      dy=2*expLambda*((MATVECT(M,y))*(MATVECT(M,phi$fun[tindex,2:ncol(phi$fun)]))) - (expLambda+mu)*y
      return(list(dy))
    }
    out <- ode(y = ini, times = timePhi, func = dKhi, parms = NULL)
    rep = out[tend,-1]
    
  }else{
    inter=min(interval,(t-s)/(phi$fun[2,1]-phi$fun[1,1]))
    do=TRUE
    while(do){
      done=TRUE
      rep = try(MagnusExpansion_ClaDS1(phi,ini,tini,tend,inter,conv=conv,timeIni=s,timeEnd=t))
      do=inherits(rep,"try-error")
      # if(! do) do=(min(rep)<0)
      if(! do & func=="Khi"){ do=(max(rep)>1)}else{do=(max(rep)>(max(1+exp(expLambda*(s-t)))*max(ini)))}# | any(rep<0))}
      inter=(inter/2)
      if(do){
        done=FALSE
        do=(inter>1e-1)
      }
    }
  }
  if(! done) rep=exp(- (expLambda+phi$mu)*(t - s))*ini
  changeRep=(! is.finite(rep) | rep<=0)
  rep[changeRep]=exp(- (expLambda[changeRep]+phi$mu)*(t - s))*ini[changeRep]#1e-310#1e-10*min(rep[is.finite(rep)])
  rep[rep==0]=1e-300
  if(func=="Psi"|func=="Zeta"){
    rep=rep/(1-phi$fun[tend,-1])
    changeRep=(! is.finite(rep) | rep<=0)
    # if(sum(is.finite(rep))>0){
    rep[changeRep]=exp(- (expLambda[changeRep]+phi$mu)*(t - s))*ini[changeRep]#1e-310#1e-10*min(rep[is.finite(rep)])
    # }else{
    #   rep[! is.finite(rep)]=1e-310
    # }
  }
  
  return(rep)
}


createLikelihood_ClaDS1 <- function(phylo, root_depth=0, relative = FALSE,nt=1000,nlambda=100,
                                    mlambda=1e-25,Mlambda=1e10,conv=1e-10){
  method="FFT"
  interval=nt
  nbtips = Ntip(phylo)
  edge = phylo$edge
  nedges= nrow(phylo$edge)
  edge.length = phylo$edge.length
  roots=c()
  div=30
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
  
  
  ll<-function(param,f,former=NULL){
    sigma=param[1]
    alpha=param[2]
    mu=param[3]
    lambda=param[-(1:3)]
    lambda2=lambda
    la=log(alpha)
    ls=log(1+sigma)
    if(relative){lambda2=relToAbs(lambda2)}
    explambda=exp(lambda2)
    Mlambda2=min(log(Mlambda),(max(lambda2)+1*max(sigma,la)))
    # Mlambda2=exp(log(2)*ceiling(log(Mlambda2)/log(2)))
    Mlambda2=exp(ls*ceiling((Mlambda2)/(ls)))
    
    mlambda2=max(log(mlambda),(min(lambda2)+1*min(-1*sigma,la)))
    # mlambda2=exp(log(2)*floor(log(mlambda2)/log(2)))
    mlambda2=exp(ls*floor((mlambda2)/(ls)))
    
    nlambda2=ceiling(max(50,min(nlambda,div*((log(Mlambda2)-log(mlambda2))/max(0.01,sigma)))))#,min(nlambda,10*((log(Mlambda2)-log(mlambda2))/abs(log(alpha))))))
    if (any(explambda>Mlambda2) | any(explambda<mlambda2)) return(list(lambda=lambda2,mu=mu,sigma=sigma,alpha=alpha,phi=NULL,PsiKhi=NULL,f=f,LL=-Inf,Pr=-Inf,LP=-Inf))
    if (any(c(mu,f,sigma,alpha) < 0)) return(list(lambda=lambda2,alpha=alpha,mu=mu,sigma=sigma,phi=NULL,PsiKhi=NULL,f=f,LL=-Inf,Pr=-Inf,LP=-Inf))
    eval=is.null(former)
    # if(!eval){eval=((abs(former$Mlambda2 - Mlambda2)>0) | (abs(former$mlambda2 - mlambda2)>0) |(abs(former$sigma - sigma)>0)|(abs(former$alpha - alpha))>0|(abs(former$mu - mu))>0|(abs(former$f - f))>0)}
    if(!eval){eval=((abs(log(former$Mlambda2) - log(Mlambda2))>0) | (abs(log(former$mlambda2) - log(mlambda2))>0) |(abs(log(former$sigma) - log(sigma))>0)|(abs(log(former$alpha) - log(alpha)))>0|(abs(former$mu - mu))>0|(abs(former$f - f))>0)}
    if(eval){
      PsiKhi=rep(0,nedges+1)
      phi=Phi_ClaDS1(sigma,alpha,mlambda2,Mlambda2,nlambda2,mu,f,tf=tf,by=tf/nt)
      for(i in 1:(length(lambda2)-1)){
        if(type[i]==1){
          m=Khi_ClaDS1(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,func = "Khi",interval=interval,conv=conv)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
        }else{
          m=Khi_ClaDS1(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                       lambda2 = lambda2[offspring[2,i]+1],func = "Ksi",nt=nt,method=method,interval=interval,conv=conv)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
        }} 
      if(root_depth==0){
        ind1=which.min(abs(lambda2[roots[1]+1]-phi$lambda))
        ind2=which.min(abs(lambda2[roots[2]+1]-phi$lambda))
        PsiKhi[1]=sum(dnorm(lambda2[roots+1],mean=lambda2[1],sd=sigma,log=TRUE))-log(max(1-phi$fun[nrow(phi$fun),ind1+1],1e-3))-log(max(1-phi$fun[nrow(phi$fun),ind2+1],1e-3))
      }else{
        m=Khi_ClaDS1(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
                     lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,interval=interval,conv=conv)
        ind=which.min(abs(lambda2[1]-phi$lambda))
        PsiKhi[1]=log(m[ind])
      }
    }else{
      change=which(abs(lambda2-former$lambda)>0)-1
      change=unique(c(change,sapply(change, function(x){if(x>0){parents[x]}else{-1}})))
      change=change[change>(-1)]
      PsiKhi=former$PsiKhi
      phi=former$phi
      for( i in change){
        if(i ==0){
          if(root_depth==0){
            ind1=which.min(abs(lambda2[roots[1]+1]-phi$lambda))
            ind2=which.min(abs(lambda2[roots[2]+1]-phi$lambda))
            PsiKhi[1]=sum(dnorm(lambda2[roots+1],mean=lambda2[1],sd=sigma,log=TRUE))-log(max(1-phi$fun[nrow(phi$fun),ind1+1],1e-3))-log(max(1-phi$fun[nrow(phi$fun),ind2+1],1e-3))
          }else{
            m=Khi_ClaDS1(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
                         lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,method=method,interval=interval,conv=conv)
            ind=which.min(abs(lambda2[1]-phi$lambda))
            PsiKhi[1]=log(m[ind])
          }
        }else{
          if(type[i]==1){
            m=Khi_ClaDS1(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,func = "Khi",interval=interval,conv=conv)
            ind=which.min(abs(lambda2[i+1]-phi$lambda))
            PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
          }else{
            m=Khi_ClaDS1(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                         lambda2 = lambda2[offspring[2,i]+1],func = "Ksi",nt=nt,method=method,interval=interval,conv=conv)
            ind=which.min(abs(lambda2[i+1]-phi$lambda))
            PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
          }
        }
      }
    }
    LL=sum(PsiKhi)
    LPr=0
    rep=list(lambda=lambda2,Mlambda2=Mlambda2,mlambda2=mlambda2,mu=mu,alpha=alpha,sigma=sigma,phi=phi,PsiKhi=PsiKhi,f=f,LL=LL,Pr=LPr,LP=LL+LPr)
    return(rep)
  }
  
  return(list(ll=ll, relToAbs=relToAbs))
  
}


#### for ClaDS2 ####


Phi_ClaDS2=function(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1,method="Higham08.b"){
  lambdaIs = seq(log(mlambda),log(Mlambda),length=nlambda+1)
  
  
  firstRow = -(lambdaIs-log(mlambda))+log(alpha)
  firstCol = lambdaIs-log(mlambda)+log(alpha)
  if(sigma==0){
    row_indices = which.min(abs(firstRow))
    firstRow = rep(0,nlambda+1)
    firstRow[row_indices] = 1
    col_indices = which.min(abs(firstCol))
    firstCol = rep(0,nlambda+1)
    firstCol[col_indices] = 1
  }else{
    firstRow = exp(-firstRow^2/(2*sigma^2))
    firstCol = exp(-firstCol^2/(2*sigma^2))
  }
  M = MPhiFFT(firstCol, firstRow, method="FFT")
  
  ini=rep((1-f),nlambda+1)
  ind=(1:length(lambdaIs))
  expLambda=exp(lambdaIs)
  
  MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
  dPhi=function(t,y,parms){
    
    dy=expLambda*(MATVECT(M, y)^2-y)+mu*expLambda*(1-y)
    return(list(dy))
  }
  
  times <- seq(from = tini, to =tf, by = by)
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL, method="daspk", atol=1e-3)
  out[,-1][out[,-1]>(1-1e-3)]=1
  out[,-1][out[,-1]<0]=0
  
  return(list(lambda=lambdaIs,expLambda=expLambda,fun=out,mu=mu,f=f,sigma=sigma,M=M,func="Phi",alpha=alpha,ind=ind))
}

#  Solve an ODE dY/dt = A(t) Y(t) by using Magnus expansion schemes
#  Order 4 with equispaced points:
#    B1 = B(Tn), B2 = B(Tn + h/2), B3 = B(Tn+h)
#    Omega(h) = (h/6)*(B1 + 4 B2 + B3) - (h*h/12) [B1,B3]
#    Y(n+1) = exp(Omega(h)) Y(n)

MagnusExpansion_ClaDS2=function(phi,ini,tini,tend,step,conv=1e-10,mask_val=-20,timeIni,timeEnd){
  expLambda=phi$expLambda
  mu=phi$mu*expLambda
  timePhi=phi$fun[,1]
  M=phi$M
  iniRep=ini
  if(tend<(tini+step)) tend=tini+step
  MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
  tStep=step*diff(timePhi)[1]
  S=seq(timeIni, max(timeIni,timeEnd-tStep), by=tStep)
  if(length(S)==1) {
    h=timeEnd-timeIni
  }else{
    h = S[2]-S[1]
  }
  D0 = expLambda+mu
  mask = - (expLambda+mu)*h < mask_val
  D3 = expLambda*MATVECT(M, phi$fun[tini+step,-1])
  for (i in S) {
    # Compute iteratively exp(Omega)*V
    
    # Magnus expansion of order 4 with equispaced points:
    #    B1 = B(Tn), B2 = B(Tn+h/2), B3 = B(Tn+h)
    #    Omega(h) = (h/6)*(B1 + 4*B2 + B3) - (h*h/12) [B1,B3]
    #    Y(n+1) = exp(Omega(h)) Y(n)
    #
    # We compute Y(n+1) = exp(Omega(h)) %*% Y(n) by an iterative process
    #    OV(k+1) = Omega %*% OV(k) / (k+1)
    #    EXPOV(k+1) = EXPOV(k) + OV(k+1)
    # Series EXPOV(k) converges towards exp(Omega(h)) %*% Y(0), and we only need
    # to explicit the product of Omega(h) by a vector.
    #
    # Here, B(t) = -diag(expLambda+mu) + (2*expLambda*diag(M %*% Phi(t))) * M
    #    D0 = diag(expLambda+mu)
    #    D1 = expLambda*diag(M %*% Phi(Tn)),
    #    D2 = expLambda*diag(M %*% Phi(Tn+h/2))
    #    D3 = expLambda*diag(M %*% Phi(Tn+h))
    #    D123 = D1+4*D2+D3 = expLambda*diag(M %*% (Phi(Tn) + 4 Phi(Tn+h/2) + Phi(Tn+h)))
    #
    # This gives
    #    Omega(h)  = - h*D0 + (h/3)*D123 * M - (h*h/12) [B1,B3]
    #
    # First two terms are straightforward, we now have to explicit the last one  [B1,B3] = B1*B3 - B3*B1
    # Recall that
    #    B1 = - D0 + 2*D1*M
    #    B3 = - D0 + 2*D3*M
    # Thus,
    #    B1 * B3 = D0^2 - 2*D0*D3*M - 2*D1*M*D0 + 4*D1*M*D3*M
    #    B3 * B1 = D0^2 - 2*D0*D1*M - 2*D3*M*D0 + 4*D3*M*D1*M
    #    [B1,B3] = - 2*D0*(D3-D1)*M + 2*(D3-D1)*M*D0 + 4*(D1*M*D3*M - D3*M*D1*M)
    #            = - 2*(D3-D1)*(D0*M - M*D0) + 4*(D1*M*D3*M - D3*M*D1*M)
    # Eventually,
    #    Omega(h) = - h*D0 + (h/3)*D123*M - (h*h/12)*[B1,B3]
    #             = - h*D0 + (h/3)*D123*M + (h*h/6)*(D3-D1)*(D0*M - M*D0) - (h*h/3)*(D1*M*D3*M - D3*M*D1*M)
    
    if(length(S)==1){
      phi123 = InterpolatedPhi(phi,timeIni) + 4 * InterpolatedPhi(phi,(timeEnd+timeIni)/2) + InterpolatedPhi(phi,timeEnd)
      D1 = D3
      D3 = expLambda*MATVECT(M, InterpolatedPhi(phi,timeEnd))
    }else{
      phi123 = InterpolatedPhi(phi,i)+ 4 * InterpolatedPhi(phi,i+tStep/2) + InterpolatedPhi(phi,i+tStep)
      D1 = D3
      D3 = expLambda*MATVECT(M,InterpolatedPhi(phi,i+tStep))
    }
    D123 = expLambda*MATVECT(M,phi123)
    
    ini[mask] = 0
    EXPOV=ini
    OV=ini
    k=1
    epsnormv = conv * sum(abs(OV))
    end=0
    while (sum(abs(OV)) > epsnormv & end < 1e3){
      end=end+1
      MOV = MATVECT(M,OV)
      OV = - h*D0*OV + (h/3)*D123*MOV + (h*h/6)*(D3-D1)*(D0*MOV - MATVECT(M,D0*OV)) - (h*h/3)*(D1*MATVECT(M,D3*MOV)-D3*MATVECT(M,D1*MOV))
      OV = OV/k
      OV[mask] = 0
      EXPOV = EXPOV + OV
      k=k+1
    }
    if(end==1e3){EXPOV[]=0}
    EXPOV[mask] = 0
    ini = EXPOV
  }
  EXPOV[mask] = exp(- (expLambda[mask]+mu[mask])*(timeEnd - timeIni))*iniRep[mask]
  return(EXPOV)
}

Khi_ClaDS2=function(phi,s,t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
                    timePhi=phi$fun[,1],nt=1000,method="Higham08.b",sparse=FALSE,interval=2,conv=1e-10){
  
  if((0.99*t)>phi$fun[nrow(phi$fun),1]) return(rep(1e-300,length(phi$expLambda)))
  if( !(func %in% c("Khi","Psi","Ksi","Zeta"))){stop("the function is not known")}
  if(any(is.na(phi$fun))){stop("Phi is NA")}
  
  tini=which.min(abs(timePhi-s))
  tend=which.min(abs(timePhi-t))
  nlambda=length(lambdas)
  expLambda=phi$expLambda
  
  if(func=="Khi" | func=="Psi"){
    
    ini=1-phi$fun[tini,-1]  
    
  }else{
    #on transforme la proba d'être dans un intervalle (de taille lambdas[2]-lambdas[1]) en densité 
    ini=expLambda*dnorm(lambda1,lambdas+log(phi$alpha),phi$sigma)*dnorm(lambda2,lambdas+log(phi$alpha),phi$sigma)
    
  }
  
  if(method == "ode"){
    MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
    dKhi=function(t,y,parms){
      tindex=which.min(abs(timePhi-t))
      dy=2*expLambda*((MATVECT(M,y))*(MATVECT(M,phi$fun[tindex,2:ncol(phi$fun)]))) - (expLambda*(1+mu))*y
      return(list(dy))
    }
    out <- ode(y = ini, times = timePhi, func = dKhi, parms = NULL)
    rep = out[tend,-1]
    
  }else{
    inter=min(interval,(t-s)/(phi$fun[2,1]-phi$fun[1,1]))
    do=TRUE
    while(do){
      done=TRUE
      rep = try(MagnusExpansion_ClaDS2(phi,ini,tini,tend,inter,conv=conv,timeIni=s,timeEnd=t))
      do=inherits(rep,"try-error")
      # if(! do) do=(min(rep)<0)
      if(! do & func=="Khi"){ do=(max(rep)>1)}else{do=(max(rep)>(max(1+exp(expLambda*(s-t)))*max(ini)))}# | any(rep<0))}
      inter=(inter/2)
      if(do){
        done=FALSE
        do=(inter>1e-1)
      }
    }
  }
  if(! done) rep=exp(- (expLambda*(1+phi$mu))*(t - s))*ini
  changeRep=(! is.finite(rep) | rep<=0)
  rep[changeRep]=exp(- (expLambda[changeRep]*(1+phi$mu))*(t - s))*ini[changeRep]#1e-310#1e-10*min(rep[is.finite(rep)])
  rep[rep==0]=1e-300
  if(func=="Psi"|func=="Zeta"){
    rep=rep/(1-phi$fun[tend,-1])
    changeRep=(! is.finite(rep) | rep<=0)
    # if(sum(is.finite(rep))>0){
    rep[changeRep]=exp(- (expLambda[changeRep]*(1+phi$mu))*(t - s))*ini[changeRep]#1e-310#1e-10*min(rep[is.finite(rep)])
    # }else{
    #   rep[! is.finite(rep)]=1e-310
    # }
  }
  
  return(rep)
}


createLikelihood_ClaDS2 <- function(phylo, root_depth=0, relative = FALSE,nt=1000,nlambda=100,
                                    mlambda=1e-25,Mlambda=1e10,conv=1e-10){
  method="FFT"
  nbtips = Ntip(phylo)
  edge = phylo$edge
  nedges= nrow(phylo$edge)
  edge.length = phylo$edge.length
  roots=c()
  interval=nt
  div=30
  
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
  
  
  ll<-function(param,f,former=NULL){
    sigma=param[1]
    alpha=param[2]
    mu=param[3]
    lambda=param[-(1:3)]
    lambda2=lambda
    la=log(alpha)
    ls=log(1+sigma)
    if(relative){lambda2=relToAbs(lambda2)}
    explambda=exp(lambda2)
    Mlambda2=min(log(Mlambda),(max(lambda2)+1*max(sigma,la)))
    # Mlambda2=exp(log(2)*ceiling(log(Mlambda2)/log(2)))
    Mlambda2=exp(ls*ceiling((Mlambda2)/(ls)))
    
    mlambda2=max(log(mlambda),(min(lambda2)+1*min(-1*sigma,la)))
    # mlambda2=exp(log(2)*floor(log(mlambda2)/log(2)))
    mlambda2=exp(ls*floor((mlambda2)/(ls)))
    
    nlambda2=ceiling(max(50,min(nlambda,div*((log(Mlambda2)-log(mlambda2))/max(0.01,sigma)))))#,min(nlambda,10*((log(Mlambda2)-log(mlambda2))/abs(log(alpha))))))
    if (any(explambda>Mlambda2) | any(explambda<mlambda2)) return(list(lambda=lambda2,mu=mu,sigma=sigma,alpha=alpha,phi=NULL,PsiKhi=NULL,f=f,LL=-Inf,Pr=-Inf,LP=-Inf))
    if (any(c(mu,f,sigma,alpha) < 0)) return(list(lambda=lambda2,alpha=alpha,mu=mu,sigma=sigma,phi=NULL,PsiKhi=NULL,f=f,LL=-Inf,Pr=-Inf,LP=-Inf))
    eval=is.null(former)
    # if(!eval){eval=((abs(former$Mlambda2 - Mlambda2)>0) | (abs(former$mlambda2 - mlambda2)>0) |(abs(former$sigma - sigma)>0)|(abs(former$alpha - alpha))>0|(abs(former$mu - mu))>0|(abs(former$f - f))>0)}
    if(!eval){eval=((abs(log(former$Mlambda2) - log(Mlambda2))>0) | (abs(log(former$mlambda2) - log(mlambda2))>0) |(abs(log(former$sigma) - log(sigma))>0)|(abs(log(former$alpha) - log(alpha)))>0|(abs(former$mu - mu))>0|(abs(former$f - f))>0)}
    if(eval){
      PsiKhi=rep(0,nedges+1)
      phi=Phi_ClaDS2(sigma,alpha,mlambda2,Mlambda2,nlambda2,mu,f,tf=tf,by=tf/nt)
      for(i in 1:(length(lambda2)-1)){
        if(type[i]==1){
          m=Khi_ClaDS2(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,func = "Khi",interval=interval,conv=conv)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
        }else{
          m=Khi_ClaDS2(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                       lambda2 = lambda2[offspring[2,i]+1],func = "Ksi",nt=nt,method=method,interval=interval,conv=conv)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
        }} 
      if(root_depth==0){
        ind1=which.min(abs(lambda2[roots[1]+1]-phi$lambda))
        ind2=which.min(abs(lambda2[roots[2]+1]-phi$lambda))
        PsiKhi[1]=sum(dnorm(lambda2[roots+1],mean=lambda2[1],sd=sigma,log=TRUE))-log(max(1-phi$fun[nrow(phi$fun),ind1+1],1e-3))-log(max(1-phi$fun[nrow(phi$fun),ind2+1],1e-3))
      }else{
        m=Khi_ClaDS2(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
                     lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,interval=interval,conv=conv)
        ind=which.min(abs(lambda2[1]-phi$lambda))
        PsiKhi[1]=log(m[ind])
      }
    }else{
      change=which(abs(lambda2-former$lambda)>0)-1
      change=unique(c(change,sapply(change, function(x){if(x>0){parents[x]}else{-1}})))
      change=change[change>(-1)]
      PsiKhi=former$PsiKhi
      phi=former$phi
      for( i in change){
        if(i ==0){
          if(root_depth==0){
            ind1=which.min(abs(lambda2[roots[1]+1]-phi$lambda))
            ind2=which.min(abs(lambda2[roots[2]+1]-phi$lambda))
            PsiKhi[1]=sum(dnorm(lambda2[roots+1],mean=lambda2[1],sd=sigma,log=TRUE))-log(max(1-phi$fun[nrow(phi$fun),ind1+1],1e-3))-log(max(1-phi$fun[nrow(phi$fun),ind2+1],1e-3))
          }else{
            m=Khi_ClaDS2(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
                         lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,method=method,interval=interval,conv=conv)
            ind=which.min(abs(lambda2[1]-phi$lambda))
            PsiKhi[1]=log(m[ind])
          }
        }else{
          if(type[i]==1){
            m=Khi_ClaDS2(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,func = "Khi",interval=interval,conv=conv)
            ind=which.min(abs(lambda2[i+1]-phi$lambda))
            PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
          }else{
            m=Khi_ClaDS2(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                         lambda2 = lambda2[offspring[2,i]+1],func = "Ksi",nt=nt,method=method,interval=interval,conv=conv)
            ind=which.min(abs(lambda2[i+1]-phi$lambda))
            PsiKhi[i+1]=InterpolatedKhi(phi,m,lambda2[i+1])
          }
        }
      }
    }
    LL=sum(PsiKhi)
    LPr=0
    rep=list(lambda=lambda2,Mlambda2=Mlambda2,mlambda2=mlambda2,mu=mu,alpha=alpha,sigma=sigma,phi=phi,PsiKhi=PsiKhi,f=f,LL=LL,Pr=LPr,LP=LL+LPr)
    return(rep)
  }
  
  return(list(ll=ll, relToAbs=relToAbs))
  
}



#### the fit_ClaDS function ####

fit_ClaDS = function(tree,sample_fraction,iterations, thin = 50, file_name = NULL, it_save = 1000,
                     model_id="ClaDS2", nCPU = 1, mcmcSampler = NULL, verbose=TRUE, ...){
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
        nCPU = nCPU,
        verbose=verbose)   
      
      if (! is.null(file_name)){
        saveRDS(mcmcSampler, file = file_name)
      }
    }
  }
  
  if((n_run * it_save)< iterations){
    mcmcSampler = add_iter_ClaDS(           
      mcmcSampler = mcmcSampler,     
      iterations = iterations - n_run *it_save,            
      thin = thin,  
      nCPU = nCPU,
      verbose=verbose)   
    
    if (! is.null(file_name)){
      saveRDS(mcmcSampler, file = file_name)
    }
  }

  return(mcmcSampler)
  
}
