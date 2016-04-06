library(deSolve)
library(fields)
library(expoRkit)
library(rexpokit)
library(Matrix)

source("birthdeath.tree.rateshift.R")

######functions#############

#Proba de ne pas avoir de descendant dans la phylogenie reconstruite 
Phi=function(sigma,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1){
  lambdaIs=seq(0,Mlambda,Mlambda/nlambda)
  if(sigma==0){
    M=diag(rep(1,nlambda+1))
  }else{
    M=t(sapply(lambdaIs,function(e){
      s=sapply(lambdaIs,function(e2){(1/(sqrt(6.28)*sigma))*exp(-(e-e2)^2/(2*sigma^2))});return(s/sum(s))
    }))}
  ini=rep((1-f),nlambda+1)
  
  dPhi=function(t,y,parms){
    dy=lambdaIs*((M%*%y)^2-y)+mu*(1-y)
    return(list(dy))
  }
  
  times <- seq(from = tini, to =tf, by = by)
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL)
  
  return(list(lambda=lambdaIs,fun=out,mu=mu,f=f,sigma=sigma,M=M,func="Phi"))
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prodMatVect_FFT=function(G, v, v0){
    w=Re(fft(G*fft(c(v,v0)),inverse=TRUE)[0:length(v)])/length(G) 
    return(w)
}

norm_vect1=function(x) sum(abs(x))

expMatBVect_FFT=function(G, v, D1, D2, D3, dt, eps){
    mask_D1 = D1 < -20/dt
    v[mask_D1] = 0
    expGv=v
    Gnv=v
    n=1
    N= length(v)
    NG = length(G)
    v0 = rep(0,N)
    D4=D3*D2
    epsnormv = eps*norm_vect1(v)
    
    while (norm_vect1(Gnv) > epsnormv){ 
      Gnv = D1*Gnv + D4*Re(fft(G*fft(c(Gnv,v0)),inverse=TRUE)[0:N])/NG
      Gnv = (dt/n)*Gnv
      Gnv[mask_D1] = 0
      expGv = expGv + Gnv
      n=n+1
      }
      expGv[mask_D1] = 2e-9
    return(expGv)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Proba de ne pas avoir de descendant dans la phylogenie reconstruite 
Phi_FFT=function(sigma,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1){

  coeff = Mlambda/nlambda
  lambdaIs  = seq(0,nlambda,length=nlambda+1)*coeff
  lambdaIs2 =-seq(nlambda+1,1,length=nlambda+1)*coeff
  lambdaIsT = c(lambdaIs,lambdaIs2)
  if(sigma==0){
    vec_density = rep(1,nlambda+1)
  }else{
    vec_density = exp(-lambdaIsT^2/(2*sigma^2))
  } 
  G=Re(fft(vec_density))
  # ------- Toeplitz matrix -------------
  M = toeplitz(exp(-lambdaIs^2/(2*sigma^2)))
  normM = 1/rowSums(M)
  M = normM * M 

  ini=rep((1-f),nlambda+1)
  vect0 = rep(0,length(ini))  
  
  dPhi=function(t,y,parms){
    My=normM*prodMatVect_FFT(G,y,vect0)
    dy=lambdaIs*(My^2-y)+mu*(1-y)
    return(list(dy))
  }

  times <- seq(from = tini, to =tf, by = by)
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL)

  return(list(lambda=lambdaIs,fun=out,mu=mu,f=f,sigma=sigma,M=M,func="Phi",G=G,normM=normM))
}

plot.Phi_FFT=function(rep,lambda,xleg=1,yleg=0,ylim=c(0,1),legend=3){
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

#calcul des fonctions Khi, Psi, Ksi, Zeta (cf document). Le calcul de ces quatres fonctions utilsent phi, donc cette fonction 
#prend en entrée la sortie de la fonction Phi (pour ne pas la recalculer dans le calcul de la vraisemblance)

Khi=function(phi,s,t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
                      timePhi=phi$fun[,1],nt=1000,method="Higham08.b",banded=0,sparse=F){

  if(t>phi$fun[nrow(phi$fun),1]) stop("t must be lower than phi$t")
  if( !(func %in% c("Khi","Psi","Ksi","Zeta"))){stop("the function is not known")}
  
  tini=which.min(abs(timePhi-s))
  tend=which.min(abs(timePhi-t))
  nlambda=length(lambdas)
  
  if(func=="Khi" | func=="Psi"){
    
    ini=1-phi$fun[tini,-1]  
    
  }else{
    #on transforme la proba d'être dans un intervalle (de taille lambdas[2]) en densité 
    ind=c(which.min(abs(lambdas-lambda1)),which.min(abs(lambdas-lambda2)))
    temp = (1-phi$fun[tini,ind[1]+1])*(1-phi$fun[tini,ind[2]+1])/(lambdas[2])^2  
    ini = temp*(M[ind[1],]*M[ind[2],]*lambdas)
    # print(ini)
  }
  
  if (substring(method,1,3) == "FFT") {
    vect0 = rep(0,length(ini)) 
    if(tini==tend){
      A = as.vector(2*lambdas*(phi$normM*prodMatVect_FFT(phi$G,phi$fun[tini,-1],vect0)))
    }else{
      SumPhi=colSums(phi$fun[tini:tend,-1])
      A = as.vector(2*lambdas*(phi$normM*prodMatVect_FFT(phi$G,SumPhi,vect0)))
      A=A/(tend-tini+1)
    }
    if(method=="FFT_none"){
      rep=rep(1, length(ini))
    }else{
      rep = expMatBVect_FFT(phi$G,ini,-(lambdas+mu), phi$normM, A,t-s,1e-20)
    }
  }else{
    if(tini==tend){
      A=as.vector(2*lambdas*(M %*% phi$fun[tini,-1]))
      B=(-(diag(lambdas+mu))+M*A)*(t-s)
    }else{
      A=as.vector(2*lambdas*as.vector(M %*% rowSums(t(phi$fun[tini:tend,-1]))))
      B=(-(diag(lambdas+mu))+ M*A /(tend-tini+1))*(t-s)}
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


library(parallel)

#fonction qui retourne le log de la fonction de vraisemblance
#les taux de naissances sont donnés dans le vecteur lambda : le premier élément est le taux de speciation initial,
#les autres sont les différences entre le taux de specation d'une branche et celui de la branche parente.
#La fonction relToAbs transforme ce vecteur en un vecteur dont le premier élément est le taux initial de speciation,
#et le ieme élement (i>2) est le taux de speciation de la (i-1)eme branche de l'arbre

createLikelihood_death <- function(phylo, root_depth=0, relative = F,nt=1000,nlambda=100,n=20, 
                                            nCPU=1, method="Higham08.b",banded=0,sparse=F){
  
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
    ll<-function(lambda, sigma, mu,f){
      lambda2=lambda
      lambda2=relToAbs(lambda2)
      
      if (any(lambda2 <= 0)) return(-Inf)
      if (any(c(mu,f,sigma) < 0)) return(-Inf)
      M=max(lambda2)*(1+n*sigma)
  
      phi=Phi_FFT(sigma,M,nlambda,mu,f,tf=tf,by=tf/nt)
  
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
    ll<-function(lambda, sigma, mu,f){
      lambda2=lambda
      lambda2=relToAbs(lambda2)
      
      if (any(lambda2 <= 0)) return(-Inf)
      if (any(c(mu,f,sigma) < 0)) return(-Inf)
      M=max(lambda2)*(1+n*sigma)
  
      phi=Phi(sigma,M,nlambda,mu,f,tf=tf,by=tf/nt)
  
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

if(T){
  seed=92146169
  set.seed(seed)
  d=0
  a=birthdeath.tree.rateshift(1,0.1,d,sigma=0.03,taxa.stop = 100,condition="taxa",new_lamb_law = "normal",mu_min = d,mu_max=d,prune.extinct = T)
  tree=a$tree
  ntips=tree$Nnode+1
  nedges=2*tree$Nnode
  size_nt = 800
  
  true.rate=a$lamb$par[a$rates]
  lambda=abs.to.rel(tree,true.rate,0.1)
  # lambda=c(1,runif(nedges,0,0.1))
  sigma=0.1
  mu=0.02
  f=1
  nCPU=1

  t2 <- Sys.time()
  name_method="FFT" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nt,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,mu,f)=[",format(c(sigma,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
 
  t2 <- Sys.time()
  name_method="FFT_none" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nt,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,mu,f)=[",format(c(sigma,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
 
  t2 <- Sys.time()
  name_method="expoRkit" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nt,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,mu,f)=[",format(c(sigma,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
 
  t2 <- Sys.time()
  name_method="none" 
  fun=createLikelihood_death(tree,relative = T,nlambda = size_nt,nt = size_nt,method=name_method , nCPU=nCPU) 
  res=fun(lambda,sigma,mu,f)
  t3 <- Sys.time()
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("(sigma,mu,f)=[",format(c(sigma,mu,f)),"] fun:", format(res))
  cat(" method=", format(name_method)," N=",format(size_nt)," CPU time:", format(t3-t2), "\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
 
}
  
