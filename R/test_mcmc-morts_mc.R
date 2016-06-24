library(TESS)
library(coda)
source("mcmcDE_mc_gibbs.R")
setwd("/Users/maliet/Dropbox/ENS/Taux/pourCluster/")
source("checkConvergence.R")

d=0.05

  seed=78877
  
  set.seed(seed)
  a=birthdeath.tree.rateshift(1,0.1,d,sigma=0.06,taxa.stop = 20,condition="taxa",new_lamb_law = "normal",mu_min = d,mu_max=d,prune.extinct = T)
  tree=a$tree
  true.rate=a$lamb$par[a$rates]
  true.rate.div=a$lamb$par[a$rates]-a$mu$par[a$rates]
  ntips=tree$Nnode+1
  nedges=2*tree$Nnode
  c=(1:nedges)[tree$edge[,1]==(ntips+1)]
  plot.with.rate(tree,true.rate,lwd=4)
  # title(main=i)
  
  
  times=as.numeric(branching.times(tree))
  target=function(x) -tess.likelihood(times,x,0)
  mini=optimize(target,c(0,10))$minimum
  
  start=list()
  N=1
  start=lapply(1:N,function(i){c(sigma=1e-4,mu=0.1,lambda_0=mini,runif(nedges,0,0))})
  # start=lapply(1:N,function(i){c(sigma=0.0001,mu=0,lambda_0=mini,runif(nedges,0,0))})
  
  trans=function(x){c(x[2:(c[2]+1)],-1*x[c[1]+2],x[(c[2]+2):(nedges+1)])}
  
  likelihood_relative=createLikelihood_death_gibbs(tree,relative = T,nt = 1000,M=2,banded = 0,method = "FFT")
  target <- function(x,former=NULL) {if(x[1]<1e-4){-Inf}else{likelihood_relative(x[-(1:2)],x[1],x[2] ,1,former=former)}}
  a=target(start[[1]])
  g=5
  testGenerator <- proposalGeneratorFactoryDE_gibbs(decreasing.var=c(rep(1e-5,3),rep(1e-5,nedges)),alpha=log(0.1)/100,
                                                    var=c(rep(1e-7,3),rep(1e-5,nedges)),
                                                    proba.gibbs=c(rep(1,g-1),1,rep(0,nedges+3-g)),
                                                    burn = 1/2,N.sigma=0,p=0.1)
  sampler=mcmcSamplerDE_gibbs(target,startvalue = start,proposalGenerator = testGenerator,Nchain = N,consoleupdates = 1)
  
  sampler=getSamplesDE_gibbs(mcmcSampler = sampler,iterations = 10,thin = 10,nCPU = 1)

  plot(sampler$codaChain)
  