
source("VraisemblanceAlphaPsi.R")
source("birthdeath.tree.rateshift.R")


check.Psi=function(t,sigma,alpha,mu,lambda_0,f,Nsim,mlambda=0.0001,Mlambda=100,nlambda=1000,nt=100){
  phi=Phi(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=t,by=t/nt)
  psi=Khi(phi,s=0,t=phi$fun[nrow(phi$fun),1],func="Psi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
          timePhi=phi$fun[,1],nt=1000,method="expoRkit")
  khi=Khi(phi,s=0,t=phi$fun[nrow(phi$fun),1],func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
          timePhi=phi$fun[,1],nt=1000,method="expoRkit")
  
  simulatePsi=simulateKhi=0
  n=0
  i=0
  j=0
  while (i < Nsim){
    cat("\r",i,"\r")
    tree=birthdeath.tree.rateshift(theta=1,lamb_par=lambda_0,mu=mu,lamb_shift=alpha,sigma=sigma,condition="time",new_lamb_law = "lognormal*shift",mu_min = mu,mu_max=mu,prune.extinct = T,time.stop = t,taxa.stop =500)
    if(!is.null(tree)){
      if (length(tree)==1){ntip=1}else{ntip=tree$tree$Nnode+1}
      if(ntip==500){
        simulatePsi=simulatePsi+0
        simulateKhi=simulateKhi+0
        n=n+1
      }else{
        simulatePsi=simulatePsi+((1-f)^(ntip-1)*f*ntip)
        simulateKhi=simulateKhi+((1-f)^(ntip-1)*f*ntip)
        n=n+(1-(1-f)^ntip)
      }
      
      i=i+1
      # n=n+1
      
    }else{
      n=n+0
    }
    j=j+1
  }
  simulatedPsi=simulatePsi/n
  simulatedKhi=simulateKhi/j
  simulatedPhi=1-(n/j)
  lind=which.min(abs(phi$expLambda-lambda_0))
  computedPsi=psi[lind]
  computedKhi=khi[lind]
  timePhi=phi$fun[,1]
  tind=which.min(abs(timePhi-t))
  computedPhi=phi$fun[tind,lind+1]
  return(list(simulatedKhi=simulatedKhi,simulatedPsi=simulatedPsi,computedKhi=computedKhi,computedPsi=computedPsi,
              simulatedPhi=simulatedPhi,computedPhi=computedPhi))
}


dataFramePsi=data.frame()
# load("pbAlpha.Rdata")
id=1
for(t in exp(seq(-2,6,length.out = 20))[1:16]){
  print(t)
  set.seed(floor(t))
  sigma=0.2
  alpha=0.7
  mu=0.1
  lambda_0=0.2
  f=1
  nlambda=200
  mlambda=exp(log(lambda_0) - 5*sigma)
  Mlambda=exp(log(lambda_0) + 5*sigma)
  nt=200
  
  cP=try(check.Psi(t,sigma,alpha,mu,lambda_0,f,Nsim=10000,nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda))
  
  if(! inherits(cP,"try-error")){
    if(nrow(dataFramePsi)==0){
      dataFramePsi=data.frame(t=t,sigma=sigma,alpha=alpha,mu=mu,lambda_0=lambda_0,f=f,
                              computedPsi=cP$computedPsi,simulatedPsi=cP$simulatedPsi,
                              computedPhi=cP$computedPhi,simulatedPhi=cP$simulatedPhi,
                              computedKhi=cP$computedKhi,simulatedKhi=cP$simulatedKhi,id=id,
                              nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda)
    }else{
      dataFramePsi=rbind(dataFramePsi,data.frame(t=t,sigma=sigma,alpha=alpha,mu=mu,lambda_0=lambda_0,f=f,
                                                 computedPsi=cP$computedPsi,simulatedPsi=cP$simulatedPsi,
                                                 computedPhi=cP$computedPhi,simulatedPhi=cP$simulatedPhi,
                                                 computedKhi=cP$computedKhi,simulatedKhi=cP$simulatedKhi,id=id,
                                                 nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda))
    }}
  print(dataFramePsi[nrow(dataFramePsi),])
}

func=function(x){(x)}
par(mfrow=c(2,2))
plot(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedKhi"],dataFramePsi[dataFramePsi[,"id"]==id,"computedKhi"]/dataFramePsi[dataFramePsi[,"id"]==id,"simulatedKhi"],xlab = "simulated chi", ylab = "computed chi / simulated chi")
plot(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedPhi"],dataFramePsi[dataFramePsi[,"id"]==id,"computedPhi"]/dataFramePsi[dataFramePsi[,"id"]==id,"simulatedPhi"],xlab = "simulated Phi", ylab = "computed Phi / simulated Phi")
legend(30,1,col=c(1,2),legend = c("computed","simulated"),lty = 1)
plot(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"computedKhi"]),type='l',xlab="t",ylab="chi")
legend(30,0.5,col=c(1,2),legend = c("computed","simulated"),lty = 1)
lines(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedKhi"]),col="red")
plot(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"computedPhi"]),type='l',xlab="t",ylab="Phi")
lines(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedPhi"]),col="red")


save(dataFramePsi,file="pbAlpha.Rdata")

