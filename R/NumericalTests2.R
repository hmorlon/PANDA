
source("VraisemblanceAlphaPsi.R")
source("birthdeath.tree.rateshift.R")


check.Psi=function(t,sigma,alpha,mu,lambda_0,f,Nsim,mlambda=0.0001,Mlambda=100,nlambda=1000,nt=100,method="Higham08.b"){

  Magnus_method="Magnus"
  #Magnus_method="Magnus2"
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  t0 <- Sys.time()
  phi=Phi(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=t,by=t/nt,method=method)
  t1 <- Sys.time()
  cat("Phi_FFT CPU time:", format(t1-t0),"\n")
  t0 <- Sys.time()
  phi2=Phi(sigma,alpha,mlambda,Mlambda,nlambda,mu,f,tini=0,tf=t,by=t/nt,method="NoFFT")
  t1 <- Sys.time()
  cat("Phi     CPU time:", format(t1-t0),"\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(" test Phi : |Phi_FFT-Phi|          =",format(norm(phi2$fun-phi$fun)),"\n")
  cat(" test Phi : |Phi_FFT-Phi|/|Phi_FFT|=",format(norm(phi2$fun-phi$fun)/norm(phi$fun)),"\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

  t0 <- Sys.time()
  psi=Khi(phi,s=0,t=phi$fun[nrow(phi$fun),1],func="Psi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
          timePhi=phi$fun[,1],nt=1000,method="ode")
  t1 <- Sys.time()
  cat("Psi ODE        CPU time:", format(t1-t0),"\n")
  t0 <- Sys.time()
  psi2=Khi(phi,s=0,t=phi$fun[nrow(phi$fun),1],func="Psi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
          timePhi=phi$fun[,1],nt=1000,method=Magnus_method)
  t1 <- Sys.time()
  cat("Psi Magnus_FFT CPU time:", format(t1-t0),"\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(" test Psi : |Psi_ODE-Psi|          =",format(norm(psi2-psi)),"\n")
  cat(" test Psi : |Psi_ODE-Psi|/|Psi_ODE|=",format(norm(psi2-psi)/norm(psi)),"\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  t0 <- Sys.time()
  khi=Khi(phi,s=0,t=phi$fun[nrow(phi$fun),1],func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
          timePhi=phi$fun[,1],nt=1000,method="ode")
  t1 <- Sys.time()
  cat("Khi ODE        CPU time:", format(t1-t0),"\n")
  t0 <- Sys.time()
  khi2=Khi(phi,s=0,t=phi$fun[nrow(phi$fun),1],func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
          timePhi=phi$fun[,1],nt=1000,method=Magnus_method)
  t1 <- Sys.time()
  cat("Khi Magnus_FFT CPU time:", format(t1-t0),"\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat(" test Khi : |Khi_ODE-Khi|          =",format(norm(khi2-khi)),"\n")
  cat(" test Khi : |Khi_ODE-Khi|/|Khi_ODE|=",format(norm(khi2-khi)/norm(khi)),"\n")
  cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  #cat("",format((khi2-khi)),"\n")

  lind=which.min(abs(phi$expLambda-lambda_0))
  computedPsi=psi2[lind]
  computedKhi=khi2[lind]
  simulatedPsi=psi[lind]
  simulatedKhi=khi[lind]
  timePhi=phi$fun[,1]
  tind=which.min(abs(timePhi-t))
  computedPhi=phi2$fun[tind,lind+1]
  simulatedPhi=phi$fun[tind,lind+1]
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
  mlambda=exp(log(lambda_0/sqrt(alpha)) - 12*sigma)
  Mlambda=exp(log(lambda_0/sqrt(alpha)) + 12*sigma)
  nt=200
  
  cP=try(check.Psi(t,sigma,alpha,mu,lambda_0,f,Nsim=10000,nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda, method = "FFT"))
  
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
plot(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedKhi"],dataFramePsi[dataFramePsi[,"id"]==id,"computedKhi"]/dataFramePsi[dataFramePsi[,"id"]==id,"simulatedKhi"],xlab = "chi ODE", ylab = "chi MAGNUS / chi ODE")
plot(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedPhi"],dataFramePsi[dataFramePsi[,"id"]==id,"computedPhi"]/dataFramePsi[dataFramePsi[,"id"]==id,"simulatedPhi"],xlab = "Phi FFT", ylab = "Phi / Phi FFT")
legend(30,1,col=c(1,2),legend = c("computed","simulated"),lty = 1)
plot(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"computedKhi"]),type='l',xlab="t",ylab="chi")
legend(30,0.5,col=c(1,2),legend = c("MAGNUS","ODE"),lty = 1)
lines(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedKhi"]),col="red")
plot(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"computedPhi"]),type='l',xlab="t",ylab="Phi")
lines(dataFramePsi[dataFramePsi[,"id"]==id,"t"],func(dataFramePsi[dataFramePsi[,"id"]==id,"simulatedPhi"]),col="red")


save(dataFramePsi,file="pbAlpha.Rdata")

