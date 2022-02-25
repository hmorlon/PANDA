fit_t_standard <- function(phylo, data, model=c("BM","OU","EB","BMtrend"), error=NULL, nuisance=FALSE, two.regime=FALSE, method="Nelder-Mead", echo=TRUE, ...){
  pars = NULL
  model = match.arg(model)[1]
  args = list(...)
    if(is.null(args[["upper"]])) upper <- Inf else upper <- args$upper
    if(is.null(args[["lower"]])) lower <- -Inf else lower <- args$lower
    if(is.null(args[["fixedRoot"]])) fixedRoot <- TRUE else fixedRoot <- FALSE
    if(is.null(args[["pars"]]))  pars <- NULL else pars <- args$pars
    
   if(inherits(phylo,"simmap")==TRUE) phylo <- reorderSimmap(phylo, order="postorder")
        else phylo <- reorder(phylo, order="postorder")
    
   if(two.regime &&!inherits(phylo,"simmap")){stop("phylo must be a simmap object with competitive regimes")}

  nobs <- Ntip(phylo)
  
  # Parameters
  if(!is.null(names(data))) {
  	data <- data[phylo$tip.label]
  	} else {
	stop("data missing taxa names")
  	}
    
  # Error provided?
  if(is.null(error)){
      error = 0
  }
  
  # if error set to NA=> automatically estimate nuisance
  if(is.na(error) & nuisance==FALSE){
     nuisance=TRUE
     error = 0
  }
  
  # index error
  index_error<-sapply(1:nobs, function(x){ which(phylo$edge[,2]==x)})
  
  # for OUM
  if(model=="OU" && two.regime==TRUE){
      precalc<-mv.Precalc(phylo, nb.traits=1, param=list(model="OUM", root=FALSE))
  }
  
  # for OU1 et random root
  if(model=="OU" && two.regime==FALSE){
      if(!is.ultrametric(phylo)) root = TRUE else root = FALSE
      precalc<-mv.Precalc(phylo, nb.traits=1, param=list(model="OU1", root=root))
  }
  
  if(two.regime){
  mod=paste0(model,"M")
  } else {
  mod=paste0(model,"1")
  }
  
  if(two.regime){ 
	regimes = ncol(phylo$mapped.edge)
  	if(regimes>2){stop("only two regime fits currently possible")}
  	}  
  # starting values
  
 if(mod=="EBM"){
 	res<-suppressWarnings(.fit_t_EB(phylo,data,regime.map=phylo,error=error,method=method,echo=FALSE))
	eval(parse(text=paste0("results<-list(LH = ",res$LogLik,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 5, sig2 = ",res$rates[2,1],", r1_",colnames(res$rates)[1]," = ",res$rates[1,1],", r2_",colnames(res$rates)[2]," = ",res$rates[1,2],", nuisance = ",res$error,", z0 = ",res$anc,", convergence = ",res$convergence,")")))
 	
 } else { 
 if(is.null(pars)){
  sig_est <- mvLL(phylo, data, param=list(estim=TRUE))$sigma # first guess
  switch(mod,
        "BM1"={ 
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            pars <- log(c(sig2,err))
        },
        "BMtrend1"={
            trendVar <- diag(vcv.phylo(phylo))
            trend <- .trendVarfun(data, trendVar, phylo)[2]
            err <- log(0.05*sig_est)
            sig2 <- log((1 - 0.05)*sig_est)
            pars <- c(sig2,trend,err)
        },
        "BMM"={
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            pars <- log(c(rep(sig2, regimes), err)) 
        },
        "OU1"={
            times <- branching.times(phylo)
            hlife <- log(2)/(max(times)/(2))
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            pars <- log(c(sig2, hlife, err))
        },
        "OUM"={
            times <- branching.times(phylo)
            hlife <- log(2)/(max(times)/(2))
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            pars <- log(c(sig2, hlife, err))
        },
        "EB1"={
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            slope <- 0
            pars <- c(log(sig2),slope,log(err))
            }
        )
     if(nuisance==FALSE) pars = pars[-length(pars)]
    }
    
  # log-likelihood function
    llik <- function(par, phy, data, mod, error){
        switch(mod,
              "BM1"={
                  sigma2 = exp(par[1])
                  phy$edge.length <- phy$edge.length*sigma2 # scaling for BM sigma
                  
                      if(nuisance==TRUE){
                          est_err = exp(par[2])
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + est_err
                      }else{
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                      }
                  
                  # ll computation
                  llik <- mvLL(phy, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE))
              },
              "BMtrend1"={
                  sigma2 = exp(par[1])
                  trend = par[2]
                  phy$edge.length <- phy$edge.length*sigma2 # scaling for BM sigma
                  
                  if(nuisance==TRUE){
                      est_err = exp(par[3])
                      phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + est_err
                  }else{
                      phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                  }
                  
                  # residual variation from the trend
                  resid <- data - trend*trendVar
                  
                  # ll computation
                  llik <- mvLL(phy, resid, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE))
                  
              },
              "BMM"={
                 sigma2 = exp(par[1:regimes])
                 phy$edge.length <- phy$mapped.edge%*%sigma2 # scaling for BMM sigma
                  
                      if(nuisance==TRUE){
                          est_err = exp(par[regimes+1])
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + est_err
                      }else{
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                      }
                  
                  # ll computation
                  llik <- mvLL(phy, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE)) 
              },
              "OU1"={
                  sigma2 = exp(par[1])
                  alpha = exp(par[2])
                  
                  # Tree transformation
                if(fixedRoot){ # to add later?
                    
                    if(is.ultrametric(phy)){
                        
                        phy <- .phyOU(phy, alpha)
                        phy$edge.length <- phy$edge.length*sigma2
                        
                        if(nuisance==TRUE){
                            est_err = exp(par[3])
                            phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + est_err
                        }else{
                            phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                        }
                        
                        # ll computation
                        llik <- mvLL(phy, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE))
                    }else{
                        
                        V<-.Call(C_panda_covar_ou_fixed, A=precalc$C1, alpha=alpha, sigma=sigma2)
                        
                        if(nuisance==TRUE){
                            est_err = exp(par[3])
                            diag(V) <- diag(V) + error^2 + est_err
                        }else{
                            diag(V) <- diag(V) + error^2
                        }
                        
                        # design matrix
                        W <- .Call(C_panda_weights, nterm=as.integer(nobs), epochs=precalc$epochs,
                        lambda=alpha, S=1, S1=1,
                        beta=precalc$listReg, root=as.integer(0))
                        
                        # ll computation
                        llik <- mvLL(V, data, method="rpf", param=list(D=W))
                    }
                    
                  }else{
                  
                   V<-.Call(C_panda_covar_ou_random, A=precalc$C1, alpha=alpha, sigma=sigma2)
                    
                   if(nuisance==TRUE){
                     est_err = exp(par[3])
                     diag(V) <- diag(V) + error^2 + est_err
                   }else{
                     diag(V) <- diag(V) + error^2   
                   }
                  
                  # design matrix
                  W <- .Call(C_panda_weights, nterm=as.integer(nobs), epochs=precalc$epochs,
                           lambda=alpha, S=1, S1=1, 
                           beta=precalc$listReg, root=as.integer(0))
                  
                  # ll computation
                  llik <- mvLL(V, data, method="rpf", param=list(D=W))
                  
                  }
              },
              "OUM"={
                  
                  sigma2 = exp(par[1])
                  alpha = exp(par[2])
                  
                  #phy <- phyOU(phy, alpha)
                  #phy$edge.length <- phy$edge.length*sigma2 
                  #
                  #   if(!is.null(error) & nuisance==TRUE){
                  #         phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + nuisance
                  #    }
                  
                  # covariance matrix
                  if(fixedRoot) V<-.Call(C_panda_covar_ou_fixed, A=precalc$C1, alpha=alpha, sigma=sigma2)
                   else  V<-.Call(C_panda_covar_ou_random, A=precalc$C1, alpha=alpha, sigma=sigma2)
                  
                   if(nuisance==TRUE){
                     est_err = exp(par[3])
                     diag(V) <- diag(V) + error^2 + est_err
                   }else{
                     diag(V) <- diag(V) + error^2   
                   }
                  
                  # design matrix
                  W <- .Call(C_panda_weights, nterm=as.integer(nobs), epochs=precalc$epochs,
                           lambda=alpha, S=1, S1=1, 
                           beta=precalc$listReg, root=as.integer(0))
                  
                  # ll computation
                  llik <- mvLL(V, data, method="rpf", param=list(D=W))
              },
              "EB1"={
              	
              	sigma2=exp(par[1])
              	rate=-abs(par[2])
              	
                #phy_temp = geiger::rescale(phy,model="EB",a=rate,sigsq=sigma2) # to limit number of packages dependencies
                phy_temp = transform_EB(phy, beta=rate, sigmasq=sigma2)
                
              	if(nuisance==TRUE){
                          est_err = exp(par[3])
                          phy_temp$edge.length[index_error]<-phy_temp$edge.length[index_error]+ error^2 + est_err
                      }else{
                          phy_temp$edge.length[index_error]<-phy_temp$edge.length[index_error]+ error^2
                      }
                  
                  # ll computation
                  llik <- mvLL(phy_temp, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE)) 
              	
              })
        
        return(llik)
    }
    
    # optimization
    if(echo==TRUE) message("Start optimization. Please wait...")
    estimModel <- optim(pars,
                        fn = function(par) -llik(par, phy=phylo, data=data, mod=mod, error=error)$logl,
                        method=method,
                        upper=upper,
                        lower=lower,
                        control=list(maxit=2000)
                        )
    
    # ancestral states/optimums
    if(echo==TRUE) message("Done. retrieve parameters and results...")
    theta <- as.numeric(llik(estimModel$par, phy=phylo, data=data, mod=mod, error=error)$theta)
    
    # param
    if(model=="EB1"){
    	param = exp(estimModel$par)
    	param[2] = -abs(log(param[2]))
    } else if(model=="BMtrend"){
        param = exp(estimModel$par)
        param[2] = estimModel$par[2]
    } else {
        param = exp(estimModel$par)
    }
    
    if(nuisance){
    switch(mod,
           "BM1"={names(param)=c("sigma2","nuisance")},
           "BMtrend1"={names(param)=c("sigma2","trend","nuisance")},
           "BMM"={names(param)=c(rep("sigma2", regimes),"nuisance")},
           "OU1"={names(param)=c("sigma2","alpha","nuisance")},
           "OUM"={names(param)=c("sigma2","alpha","nuisance")
                names(theta)=colnames(phylo$mapped.edge)},
           "EB1"={names(param)=c("sigma2", "slope","nuisance")}) 
        }else{
    switch(mod,
           "BM1"={names(param)=c("sigma2")},
           "BMtrend1"={names(param)=c("sigma2","trend")},
           "BMM"={names(param)=rep("sigma2", regimes)},
           "OU1"={names(param)=c("sigma2","alpha")},
           "OUM"={names(param)=c("sigma2","alpha")
                names(theta)=colnames(phylo$mapped.edge)},
           "EB1"={names(param)=c("sigma2", "slope")})  
    }
    
    LL = -estimModel$value
    nparam = length(param)+length(theta)
    # AIC
    AIC = -2*LL+2*nparam
    # AIC corrected
    AICc = AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) 
    
    # return results
    
    res <- list(logl=LL, AIC=AIC, AICc=AICc, param=param, theta=theta, nb_param=nparam, opt=estimModel)
    
    if(mod=="BM1"){eval(parse(text=paste0("results<-list(LH = ",res$logl,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 3, sig2 = ",as.numeric(res$param[1]),", nuisance = ",as.numeric(res$param[2]),", z0 = ",res$theta,", convergence = ",res$opt$convergence,")")))}
     if(mod=="BMtrend1"){eval(parse(text=paste0("results<-list(LH = ",res$logl,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 4, sig2 = ",as.numeric(res$param[1]),", trend = ",as.numeric(res$param[2]),", nuisance = ",as.numeric(res$param[3]),", z0 = ",res$theta,", convergence = ",res$opt$convergence,")")))}
    if(mod=="OU1"){eval(parse(text=paste0("results<-list(LH = ",res$logl,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 4, sig2 = ",as.numeric(res$param[1]),", alpha = ",as.numeric(res$param[2]),", nuisance = ",as.numeric(res$param[3]),", z0 = ",res$theta,", convergence = ",res$opt$convergence,")")))}
    if(mod=="EB1"){eval(parse(text=paste0("results<-list(LH = ",res$logl,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 4, sig2 = ",as.numeric(res$param[1]),", r = ",as.numeric(res$param[2]),", nuisance = ",as.numeric(res$param[3]),", z0 = ",res$theta,", convergence = ",res$opt$convergence,")")))}
    if(mod=="BMM"){eval(parse(text=paste0("results<-list(LH = ",res$logl,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 4, sig2_",colnames(phylo$mapped.edge)[1]," = ",as.numeric(res$param[1]),", sig2_",colnames(phylo$mapped.edge)[2]," = ",as.numeric(res$param[2]),", nuisance = ",as.numeric(res$param[3]),", z0 = ",res$theta,", convergence = ",res$opt$convergence,")")))}
    if(mod=="OUM"){eval(parse(text=paste0("results<-list(LH = ",res$logl,", aic = ",res$AIC,", aicc = ",res$AICc,", free.parameters = 5, sig2 = ",as.numeric(res$param[1]),", alpha = ",as.numeric(res$param[2]),", nuisance = ",as.numeric(res$param[3]),", z0_",names(res$theta)[1]," = ",as.numeric(res$theta[1]),", z0_",names(res$theta)[2]," = ",as.numeric(res$theta[2]),", convergence = ",res$opt$convergence,")")))}
    
    }
    return(res)
}
