library(deSolve)
library(fields)
library(Matrix)

source("birthdeath.tree.rateshift.R")
source("MPhiAbstract.class.R")
source("MPhiFFT.class.R")
source("MPhiNoFFT.class.R")

######functions#############

#Proba de ne pas avoir de descendant dans la phylogenie reconstruite 
Phi=function(sigma,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1,method="Higham08.b"){
  lambdaIs=seq(0,Mlambda,length=nlambda+1)
  if (substring(method,1,3) == "FFT") {
    if(sigma==0){
      vec_density = rep(1,nlambda+1)
    }else{
      vec_density = exp(-lambdaIs^2/(2*sigma^2))
    }
    M = MPhiFFT(vec_density, vec_density, method=method)
  } else {
    if(sigma==0){
      B=diag(rep(1,nlambda+1))
    }else{
      B=t(sapply(lambdaIs,function(e){
        s=sapply(lambdaIs,function(e2){exp(-(e-e2)^2/(2*sigma^2))});return(s/sum(s))
      }))
    }
    M = MPhiNoFFT(B, method=method)
  }
  ini=rep((1-f),nlambda+1)
  
  dPhi=function(t,y,parms){
    dy=lambdaIs*(applyV(M, y)^2-y)+mu*(1-y)
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


#calcul des fonctions Khi, Psi, Ksi, Zeta (cf document). Le calcul de ces quatres fonctions utilsent phi, donc cette fonction
#prend en entrée la sortie de la fonction Phi (pour ne pas la recalculer dans le calcul de la vraisemblance)

Khi=function(phi,s,t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
                      timePhi=phi$fun[,1],nt=1000,banded=0,sparse=F){

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
  
  if(tini==tend){
    A = (t-s)*as.vector(2*lambdas*applyV(M,phi$fun[tini,-1]))
  }else{
    SumPhi=colSums(phi$fun[tini:tend,-1])
    A = as.vector(2*lambdas*applyV(M,SumPhi))
    A=A*(t-s)/(tend-tini+1)
  }
  if(length(grep("none", M@method))>0) {
    rep=rep(1, length(ini))
  }else{
    rep = expV(M,-(lambdas+mu)*(t-s), A, ini)
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
  
  
  
  ll<-function(lambda, sigma, mu,f){
    lambda2=lambda
    lambda2=relToAbs(lambda2)
    
    if (any(lambda2 <= 0)) return(-Inf)
    if (any(c(mu,f,sigma) < 0)) return(-Inf)
    M=max(lambda2)*(1+n*sigma)
  
    phi=Phi(sigma,M,nlambda,mu,f,tf=tf,by=tf/nt,method=method)
  
    logLik=mclapply(1:(length(lambda2)-1),function(i){
      if(type[i]==1){
        m=Khi(phi,0,nodeprof[phylo$edge[i,1]],nt=nt,banded = banded,sparse=sparse)
        ind=which.min(abs(lambda2[i+1]-phi$lambda))
        return(log(m[ind]))
      }else{
        m=Khi(phi,nodeprof[phylo$edge[i,2]],nodeprof[phylo$edge[i,1]],lambda1 = lambda2[offspring[1,i]+1],
                       lambda2 = lambda2[offspring[2,i]+1],func = "Zeta",nt=nt,banded = banded,sparse=sparse)
        ind=which.min(abs(lambda2[i+1]-phi$lambda))
        return(log(m[ind]))}},mc.cores = nCPU)
    # print(logLik)
    logLik=sum(sapply(logLik,function(x){x}))
    if(root_depth==0){
      logLik=logLik+sum(dnorm(lambda[roots+1],mean=0,sd=sigma,log=T))
    }else{
      m=Khi(phi,root_depth+nodeprof[nbtips+1],nodeprof[nbtips+1],lambda1 = lambda2[roots[1]+1],
                     lambda2 = lambda2[roots[2]+1],func = "Zeta",nt=nt,banded = banded,sparse=sparse)
      ind=which.min(abs(lambda2[1]-phi$lambda))
      logLik=logLik+log(m$fun[nrow(m$fun),ind])
    }
  
    return(logLik)
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
  
