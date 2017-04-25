library(deSolve)
library(fields)
library(expoRkit)
library(rexpokit)
library(Matrix)
library(signal)
library(parallel)
source("birthdeath.tree.rateshift.R")
library(pracma)

######functions#############


Toeplitz2=function (x,...) 
{
  if (!is.vector(x)) 
    stop("'x' is not a vector")
  if (!missing(...)) {
    na <- length(list(...))
    warning(sprintf(ngettext(na, "extra argument %s will be disregarded", 
                             "extra arguments %s will be disregarded"), paste(sQuote(names(list(...))), 
                                                                              collapse = ", ")), domain = NA)
  }
  n <- length(x)
  n=ceiling(n/2)
  A <- matrix(raw(), n, n)
  matrix(x[col(A) - row(A) + n], n, n)
}


Phi=function(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1){
  lambdaIs=seq(log(mlambda),log(Mlambda),(log(Mlambda)-log(mlambda))/nlambda)
  if(sigma==0){
    M=diag(rep(1,nlambda+1))
  }else{
    M=t(sapply(lambdaIs,function(e){
      s=sapply(lambdaIs,function(e2){(1/(sqrt(6.28)*sigma))*exp(-(e-e2+log(alpha))^2/(2*sigma^2))});return(s/sum(s))
    }))}
  ini=rep((1-f),nlambda+1)
  #   ini[1]=0
  #   ini[length(ini)]=0
  ind=(1:length(lambdaIs))
  expLambda=exp(lambdaIs)
  
  dPhi=function(t,y,parms){
    
    dy=expLambda*((M%*%y)^2-y)+mu*(1-y)
    return(list(dy))
  }
  
  times <- seq(from = tini, to =tf, by = by)
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL)
  
  return(list(lambda=lambdaIs,expLambda=expLambda,fun=out,mu=mu,f=f,sigma=sigma,M=M,func="Phi",alpha=alpha,ind=ind))
}

prodMatVect_FFT=function(G, v, v0){
  w=Re(ifft(G*fft(c(v,v0)))[0:length(v)])
  return(w)
}

expMatBVect_FFT=function(G, v, D1, D2, D3, dt, eps){
  mask_D1 = D1 < -20/dt
  # mask_D1 = D1
  v[mask_D1] = 0
  expGv=v
  Gnv=v
  n=1
  N= length(v)
  v0 = rep(0,N)
  D4=D3*D2
  
  while (norm(Gnv) > eps*norm(v)){ 
    Gnv = D1*Gnv + D4*Re(ifft(G*fft(c(Gnv,v0)))[0:N])
    Gnv = (dt/n)*Gnv
    Gnv[mask_D1] = 0
    expGv = expGv + Gnv
    n=n+1
  }
  expGv[mask_D1] = 2e-9
  return(expGv)
}


Phi_FFT=function(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1,method ="lsoda"){
  
  coeff = (log(Mlambda)-log(mlambda))/nlambda
  lambdaIs  = seq(0,nlambda,length=nlambda+1)*coeff
  lambdaIs2 =-seq(nlambda,1,length=nlambda)*coeff
  lambdaIsT = c(lambdaIs+log(alpha),10,lambdaIs2+log(alpha))
  if(sigma==0){
    vec_density = rep(0,2*nlambda+2)
    vec_density[which.min(abs(lambdaIsT))]=1
  }else{
    vec_density = exp(-lambdaIsT^2/(2*sigma^2))
  } 
  
  
  G=(fft(vec_density))
  # ------- Toeplitz matrix -------------
  M = Toeplitz(vec_density[1:(nlambda+1)],c(vec_density[1],vec_density[(2*nlambda+2):(nlambda+3)]))
  # M = Toeplitz2(c(vec_density[(nlambda+1):1],vec_density[(2*nlambda+2):(nlambda+3)]))
  
  invNorm=rowSums(M)
  normM = 1/rowSums(M) #; normM[normM>1]=1
  # normM = min(normM)
  M = normM * M 
  
  ini=rep((1-f),nlambda+1)
  vect0 = rep(0,length(ini))  
  expLambda=exp(lambdaIs+log(mlambda))
  dPhi=function(t,y,parms){
    A=(normM*prodMatVect_FFT(G,y,vect0))
    # if(max(A)>1) debug(dPhi)
    dy=expLambda*(A^2-y)+mu*(1-y)
    return(list(dy))
  }
  
  times <- seq(from = tini, to =tf, by = by)
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL)
  
  return(list(lambda=lambdaIs+log(mlambda),expLambda=expLambda,fun=out,mu=mu,f=f,sigma=sigma,alpha=alpha,M=M,func="Phi",G=G,normM=normM))
}


plot.Phi=function(rep,lambda,xleg=1,yleg=0,ylim=c(0,1),legend=3){
  ind=sapply(lambda, function(x){which.min(abs(rep$lambda-x))})
  colors=rainbow(length(lambda))
  col=1
  plot(rep$fun[,1],rep$fun[,ind[1]+1],ylim = ylim,xlab = "time", ylab=rep$func,type='l',col=colors[col])
  if(length(lambda)>1){
    for(i in 2:length(lambda)){
      col=col+1
      lines(rep$fun[,1],rep$fun[,ind[i]+1],type='l',col=colors[col])
    }
  }
  if(legend==1){image.plot(z = c(min(lambda),max(lambda)),col = colors, horizontal=F,legend.only = T)
  }else if (legend==2){legend(xleg,y=yleg, legend=lambda,col=colors, lty=1, cex=0.8)}
}

Khi=function(phi,s,t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
             timePhi=phi$fun[,1],nt=1000,method="Higham08.b",banded=0,sparse=F){
  
  if(t>phi$fun[nrow(phi$fun),1]) stop("t must be lower than phi$t")
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
    if (T){#(substring(method,1,3) == "FFT") {
      indL=c(which.min(abs(lambdas-lambda1)),which.min(abs(lambdas-lambda2)))
      temp = (1-phi$fun[tini,indL[1]+1])*(1-phi$fun[tini,indL[2]+1])/(lambdas[2]-lambdas[1])^2  
      ini = temp*(M[indL[1],]*M[indL[2],]*expLambda)
    }else{
      indL=c(which.min(abs(lambdas-lambda1)),which.min(abs(lambdas-lambda2)))
      i1=max(min(length(expLambda),indL[1]+floor(log(alpha)/(lambdas[2]-lambdas[1]))),1)
      i2=max(min(length(expLambda),indL[2]+floor(log(alpha)/(lambdas[2]-lambdas[1]))),1)
      temp = (1-phi$fun[tini,indL[1]+1])*(1-phi$fun[tini,indL[2]+1])/(lambdas[2]-lambdas[1])^2  
      ini = temp*(M[i1,]*M[i2,]*expLambda)
    }
  }
  
  if (substring(method,1,3) == "FFT") {
    vect0 = rep(0,length(ini)) 
    if(tini==tend){
      A = as.vector(2*expLambda*(phi$normM*prodMatVect_FFT(phi$G,phi$fun[tini,-1],vect0)))
    }else{
      SumPhi=colSums(phi$fun[tini:tend,-1])
      A = as.vector(2*expLambda*(phi$normM*prodMatVect_FFT(phi$G,SumPhi,vect0)))
#       A=A/(tend-tini+1)
      A=A/(tend-tini)
    }
    if(method=="FFT_none"){
      rep=rep(1, length(ini))
    }else{
      rep = expMatBVect_FFT(phi$G,ini,-(expLambda+mu), phi$normM, A,t-s,1e-12)
    }
  }else{
    if(tini==tend){
      A=as.vector(2*expLambda*(M %*% phi$fun[tini,phi$ind+1]))
      B=(-(diag(expLambda+mu))+M*A)*(t-s)
    }else{
      A=as.vector(2*expLambda*as.vector(M %*% rowSums(t(phi$fun[tini:tend,phi$ind+1]))))
      B=(-(diag(expLambda+mu))+ M*A /(tend-tini+1))*(t-s)}
    if(banded>0){
      Bdm=min(abs(diag(B)))
      indices=abs(B)<(Bdm/banded)
      B[indices]=rep(0,sum(indices))
    }
    if(sparse){B=Matrix(B,sparse=T)}
    if(method=="Sidje98"){
      rep=expAtv(B, ini)$eAtv
    }else if(method=="expoRkit"){
      rep=expv(x=B,v=ini,t=1)
    }else if(method=="Rexpv"){
      rep=expoRkit:::Rexpv(B@x,B@i+1,B@p+1,length(ini),ini)
    }else if(method=="none"){
      rep=rep(1, length(ini))
    }else{
      rep=expm(B,method = method) %*% ini
    }
  }
  
  if(func=="Psi"|func=="Zeta"){
    rep=rep/(1-phi$fun[tend,-1])
  }
  
  return(rep)
}


createLikelihood_death <- function(phylo, root_depth=0, relative = F,nt=1000,nlambda=100,n=20, 
                                   nCPU=1, method="Higham08.b",banded=0,sparse=F,mlambda=1e-5,Mlambda=1e5){
  
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
  
  
  
  if (substring(method,1,3) == "FFT") {
    ll<-function(lambda, sigma, alpha, mu,f){
      lambda2=lambda
      lambda2=relToAbs(lambda2)
      
      # if (any(lambda2 <= 0)) return(-Inf)
      if (any(c(mu,f, alpha, sigma) < 0)) return(-Inf)
      
      phi=Phi_FFT(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tf=tf,by=tf/nt)
      
      logLik=mclapply(1:(length(lambda2)-1),function(i){
        if(type[i]==1){
          m=Khi(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,method=method,banded = banded,sparse=sparse)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          return(log(m[ind]))
        }else{
          m=Khi(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                lambda2 = lambda2[offspring[2,i]+1],func = "Zeta",nt=nt,method=method,banded = banded,sparse=sparse)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          return(log(m[ind]))}},mc.cores = nCPU)
      # print(logLik)
      logLik=sum(sapply(logLik,function(x){x}))
      if(root_depth==0){
        logLik=logLik+sum(dnorm(lambda[roots+1],mean=0,sd=sigma,log=T))
      }else{
        m=Khi(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
              lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,method = method,banded = banded,sparse=sparse)
        ind=which.min(abs(lambda2[1]-phi$lambda))
        logLik=logLik+log(m$fun[nrow(m$fun),ind])
      }
      
      return(logLik)
    }
  }else{
    ll<-function(lambda, sigma, alpha, mu,f){
      lambda2=lambda
      lambda2=relToAbs(lambda2)
      
      # if (any(lambda2 <= 0)) return(-Inf)
      if (any(c(mu,f,sigma) < 0)) return(-Inf)
      
      phi=Phi(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tf=tf,by=tf/nt)
      
      logLik=mclapply(1:(length(lambda2)-1),function(i){
        if(type[i]==1){
          m=Khi(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,func = "Psi",method=method,banded = banded,sparse=sparse)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          return(log(m[ind]))
        }else{
          m=Khi(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                lambda2 = lambda2[offspring[2,i]+1],func = "Zeta",nt=nt,method=method,banded = banded,sparse=sparse)
          ind=which.min(abs(lambda2[i+1]-phi$lambda))
          return(log(m[ind]))}},mc.cores = nCPU)
      # print(logLik)
      logLik=sum(sapply(logLik,function(x){x}))
      if(root_depth==0){
        logLik=logLik+sum(dnorm(lambda[roots+1],mean=0,sd=sigma,log=T))
      }else{
        m=Khi(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
              lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,method = method,banded = banded,sparse=sparse)
        ind=which.min(abs(lambda2[1]-phi$lambda))
        logLik=logLik+log(m$fun[nrow(m$fun),ind])
      }
      
      return(logLik)
    }
  }
  return(ll)
  
}

abs.to.rel=function(tree,true.rate,lambda_0){
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

######exemples########

if(F){
  seed=92146169
  set.seed(seed)
  d=0
  a=birthdeath.tree.rateshift(1,0.1,d,sigma=0.03,taxa.stop = 100,condition="taxa",new_lamb_law = "normal",mu_min = d,mu_max=d,prune.extinct = T)
  tree=a$tree
  ntips=tree$Nnode+1
  nedges=2*tree$Nnode
  size_nt = 800
  size_nlambda = 1000
  
  true.rate=a$lamb$par[a$rates]
  lambda=abs.to.rel(tree,log(true.rate),0.1)
  # lambda=c(1,runif(nedges,0,0.1))
  sigma=0.1
  mu=0.
  alpha=0.5
  f=1
  nCPU=1
  
  t2 <- Sys.time()
  name_method="FFT" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nlambda,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,alpha,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,alpha,mu,f)=[",format(c(sigma,alpha,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  t2 <- Sys.time()
  name_method="FFT_none" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nlambda,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,alpha,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,alpha,mu,f)=[",format(c(sigma,alpha,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  t2 <- Sys.time()
  name_method="expoRkit" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nlambda,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,alpha,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,alpha,mu,f)=[",format(c(sigma,alpha,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  t2 <- Sys.time()
  name_method="none" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nlambda,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,alpha,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,alpha,mu,f)=[",format(c(sigma,alpha,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  likelihood_relative=createLikelihood_relative_lognormal(tree)$ll
}



