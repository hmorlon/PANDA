library(deSolve)
library(fields)
library(signal)
library(parallel)
library(pracma)

source("birthdeath.tree.rateshift.R")
source("MPhiAbstract.class.R")
source("MPhiFFT.class.R")
source("MPhiNoFFT.class.R")

######functions#############


Phi=function(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=100,by=0.1,method="Higham08.b"){
  lambdaIs = seq(log(mlambda),log(Mlambda),length=nlambda+1)

  if (substring(method,1,3) == "FFT") {
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
    M = MPhiFFT(firstCol, firstRow, method=method)
  } else {
    if(sigma==0){
      B=diag(rep(1,nlambda+1))
    }else{
      B=t(sapply(lambdaIs,function(e){
        s=sapply(lambdaIs,function(e2){exp(-(e-e2+log(alpha))^2/(2*sigma^2))});return(s/sum(s))
      }))
    }
    M = MPhiNoFFT(B, method=method)
  }

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
  out   <- ode(y = ini, times = times, func = dPhi, parms = NULL)
  
  return(list(lambda=lambdaIs,expLambda=expLambda,fun=out,mu=mu,f=f,sigma=sigma,M=M,func="Phi",alpha=alpha,ind=ind))
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

#  Solve an ODE dY/dt = A(t) Y(t) by using Magnus expansion schemes
#  Order 4 with equispaced points:
#    B1 = B(Tn), B2 = B(Tn + h/2), B3 = B(Tn+h)
#    Omega(h) = (h/6)*(B1 + 4 B2 + B3) - (h*h/12) [B1,B3]
#    Y(n+1) = exp(Omega(h)) Y(n)

MagnusExpansion=function(phi,ini,tini,tend,method="order4_eq"){
  expLambda=phi$expLambda
  mu=phi$mu
  timePhi=phi$fun[,1]
  M=phi$M
  if (method=="order4_eq") {
    h = timePhi[tini+2] - timePhi[tini]
    for (i in seq(tini, tend-2, by=2)) {
      rep = expOmegaV(h,phi,i,i+2,ini)
      ini = rep 
    }
  } else if (method=="order4_eq2") {
    step = 6
    h = timePhi[tini+step] - timePhi[tini]
    norm = 0
    last = tini
    i = tini
    while (i <= tend - step) {
      phi123 = phi$fun[i,-1] + 4 * phi$fun[i+step/2,-1] + phi$fun[i+step,-1]
      D1 = expLambda*applyV(M, phi$fun[i,-1])
      D3 = expLambda*applyV(M, phi$fun[i+step,-1])
      D0 = expLambda*applyV(M,phi123)
      Mini = applyV(M,ini)
      OV = - h*(expLambda+mu)*ini + (h/3)*( D0*Mini - h*( D1*applyV(M,D3*Mini)-D3*applyV(M,D1*Mini) ) )
      OV = OV + (h*h/6)*( (D3-D1)*applyV(M,(expLambda+mu)*ini) - (expLambda+mu)* (D3-D1)*Mini )
      norm = norm + h * sum(abs(OV)) / sum(abs(ini))
      if (norm > pi) {
        if ((i - last) %% 2 == 0) {
          i = i - 2
        } else {
          i = i - 1
        }
        if (i < last + 2) {
          i = last + 2
        }
        h = timePhi[i] - timePhi[last]
        rep = expOmegaV(h,phi,last,i,ini)
        ini = rep
        norm = 0
        last = i
      } else {
        i = i + step
      }
    }
    if (last < tend-1) {
        h = timePhi[tend] - timePhi[last]
        rep = expOmegaV(h,phi,last,tend,ini)
        ini = rep
    }
  }
  return(rep)
}

expOmegaV=function(h,phi,i1,i2,x){
  # Compute iteratively exp(Omega)*V

  expLambda=phi$expLambda
  mu=phi$mu
  M=phi$M
  
  phi123 = phi$fun[i1,-1] + 4 * phi$fun[(i1+i2)/2,-1] + phi$fun[i2,-1]
  D1 = expLambda*applyV(M, phi$fun[i1,-1])
  D3 = expLambda*applyV(M, phi$fun[i2,-1])
  D0 = expLambda*applyV(M,phi123)

  mask = - (expLambda+mu) < -20
  x[mask] = 0
  EXPOV=x
  OV=x
  i=1
  epsnormv = 1e-10 * sum(abs(OV))
  while (sum(abs(OV)) > epsnormv){
    MOV = applyV(M,OV)
    OV = - h*(expLambda+mu)*OV + (h/3)*( D0*MOV - h*( D1*applyV(M,D3*MOV)-D3*applyV(M,D1*MOV) ) )
    OV = OV + (h*h/6)*( (D3-D1)*applyV(M,(expLambda+mu)*OV) - (expLambda+mu)* (D3-D1)*MOV )
    OV = OV/i 
    OV[mask] = 0
    EXPOV = EXPOV + OV
    i=i+1
  }
  EXPOV[mask] = 2e-9
  return(EXPOV)
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
    if (T){#(substring(method,1,3) == "FFT")
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

  MATVECT <- selectMethod(applyV, c(class(M), class(ini)))
  
  if(method == "ode"){
    dKhi=function(t,y,parms){
      tindex=which.min(abs(timePhi-t))
      dy=2*expLambda*((applyV(M,y))*(applyV(M,phi$fun[tindex,2:ncol(phi$fun)]))) - (expLambda+mu)*y
      return(list(dy))
    }
    out <- ode(y = ini, times = timePhi, func = dKhi, parms = NULL)
    rep = out[tend,-1]

  }else if(method == "Magnus"){
   rep = MagnusExpansion(phi,ini,tini,tend,method="order4_eq")
  }else if(method == "Magnus2"){
   rep = MagnusExpansion(phi,ini,tini,tend,method="order4_eq2")
  }else if(method == "Magnus1"){

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
  # Here, B(t) = -diag(expLambda+mu) + (2*expLambda*(M %*% Phi(t))) * M
  #    D1 = M %*% Phi(Tn), D2 = M %*% Phi(Tn+h/2), D3 = M %*% Phi(Tn+h)
  # This gives
  #    Omega(h) %*% V
  #      = - h*(expLambda+mu)*V + (h/3)*expLambda*(D1+4*D2+D3)*applyV(M,V)
  #           - (h*h/12) [B1,B3] %*% V
  # First two terms are straightforward, we now have to explicit the last one
  #    [B1,B3] %*% V = B1 %*% B3 %*% V - B3 %*% B1 %*% V
  # Recall that
  #    B1 = -diag(expLambda+mu) + 2*expLambda*D1*M
  #    B3 = -diag(expLambda+mu) + 2*expLambda*D3*M
  # Thus,
  #    B3 %*% V = -(expLambda+mu)*V + 2*expLambda*D3*applyV(M,V)
  #    B1 %*% B3 %*% V = - B1 %*% ((expLambda+mu)*V) + 2*B1 %*% (expLambda*D3*applyV(M,V))
  #       = (expLambda+mu)^2*V - 2*expLambda*D1*applyV(M,(expLambda+mu)*V)
  #       - 2*(expLambda+mu)*expLambda*D3*applyV(M,V) + 4*expLambda*D1*applyV(M,expLambda*D3*applyV(M,V))
  # By symmetry,
  #    B3 %*% B1 %*% V
  #       = (expLambda+mu)^2*V - 2*expLambda*D3*applyV(M,(expLambda+mu)*V)
  #       - 2*(expLambda+mu)*expLambda*D1*applyV(M,V) + 4*expLambda*D3*applyV(M,expLambda*D1*applyV(M,V))
  # and eventually
  #    [B1,B3] %*% V = 2*expLambda*(D3-D1)*applyV(M,(expLambda+mu)*V) + 2*(expLambda+mu)*expLambda*(D1-D3)*applyV(M,V)
  #        + 4*expLambda*D1*applyV(M,expLambda*D3*applyV(M,V))
  #        - 4*expLambda*D3*applyV(M,expLambda*D1*applyV(M,V))
  #      = 2*expLambda*(((D3-D1)*applyV(M,(expLambda+mu)*V) - (expLambda+mu)*applyV(M,V)) +
  #          2*D1*applyV(M,expLambda*D3*applyV(M,V)) - 2*D3*applyV(M,expLambda*D1*applyV(M,V)))
 
  step = 6
  h = timePhi[tini+step] - timePhi[tini]
  norm = 0
  last = tini
  i = tini
  rep = ini
  while (i < tend - step) {
    phi123 = phi$fun[i,-1] + 4 * phi$fun[i+step/2,-1] + phi$fun[i+step,-1]
    D1 = MATVECT(M, phi$fun[i,-1])
    D3 = MATVECT(M, phi$fun[i+step,-1])
    OV = h * ( - (expLambda+mu) * rep + expLambda*applyV(M, phi123)*applyV(M,rep) / 3
      - (h/6) * expLambda * (
          (D3-D1)*(applyV(M,(expLambda+mu)*rep) - (expLambda+mu)*applyV(M,rep))
           + 2*D1*applyV(M,expLambda*D3*applyV(M,rep)) - 2*D3*applyV(M,expLambda*D1*applyV(M,rep))
      ))
    norm = norm + h * sum(abs(OV)) / sum(abs(rep))
    if (norm > pi || i >= tend - 2*step) {
      if (norm > pi) {
        if ((i - last) %% 2 == 0) {
          i = i - 2
        } else {
          i = i - 1
        }
        if (i < last + 2) {
          i = last + 2
        }
      } else {
        i = tend - step
      }
      h = timePhi[i] - timePhi[last]
      phi123 = phi$fun[i,-1] + 4 * phi$fun[i+step/2,-1] + phi$fun[i+step,-1]
      D1 = MATVECT(M, phi$fun[i,-1])
      D3 = MATVECT(M, phi$fun[i+step,-1])

      # Compute iteratively exp(Omega)*V
      mask = - (expLambda+mu) < -20
      rep[mask] = 0
      EXPOV=rep
      OV=rep
      k=1
      epsnormv = 1e-10 * sum(abs(OV))
      while (sum(abs(OV)) > epsnormv){
        OV = h * ( - (expLambda+mu) * OV + expLambda*applyV(M, phi123)*applyV(M,OV) / 3
                   - (h/6) * expLambda * (
                        (D3-D1)*(applyV(M,(expLambda+mu)*OV) - (expLambda+mu)*applyV(M,OV))
                         + 2*D1*applyV(M,expLambda*D3*applyV(M,OV)) - 2*D3*applyV(M,expLambda*D1*applyV(M,OV))
                 )) / k
        OV[mask] = 0
        EXPOV = EXPOV + OV
        k=k+1
      }
      EXPOV[mask] = 2e-9

      rep = EXPOV
      norm = 0
      last = i
    } else {
      i = i + step
    }
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
  
  ll<-function(lambda, sigma, alpha, mu,f){
    lambda2=lambda
    lambda2=relToAbs(lambda2)
    
    # if (any(lambda2 <= 0)) return(-Inf)
    if (any(c(mu,f, alpha, sigma) < 0)) return(-Inf)
    
    phi=Phi(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tf=tf,by=tf/nt,method=method)
    
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



