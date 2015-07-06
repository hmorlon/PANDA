likelihood_MC_geog<-function(phylo,data,par=c(sig2,S),geography.object) #par[1]=sig2,par[2]=sterm
{
  	if(length(par)!=2){stop("par must contain two values, one for sig2 and another for S")}
  	sig2<-abs(par[1]) 
  	sterm<-ifelse(sign(par[2])==-1,par[2],-par[2])
	V<-try(.VCV.rescale.geog(phylo,sig2,0,sterm,geography.matrix))
	if(class(V)=="try-error"){return(Inf)}
	if(any(is.na(V))){
		return(Inf)
	} else{
  	op <- getOption("show.error.messages")
  	options(show.error.messages=FALSE)
	IV=try(solve(V))
  	options(show.error.messages=op)
  if(class(IV)=="try-error"){
    IV=pseudoinverse(V) 
  	if(max(IV)==0){return(Inf)}
  }
 	data<-as.matrix(data[rownames(V)])
	I<-matrix(rep(1,length(phylo$tip.label)))
	theta<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
	D<-(I%*%theta)-data[,1]
	T<-t(D)
    log_final_lik<- -0.5*(T%*%IV%*%D)-0.5*determinant(V)$modulus-0.5*(length(phylo$tip.label)*log(2*pi))
    if(is.na(log_final_lik) | is.infinite(log_final_lik)){log_final_lik=-1000000} #minus here because after we change the sign for minimizing the loglik
    return(-log_final_lik) 
    }
}


