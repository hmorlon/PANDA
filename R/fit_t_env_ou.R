################################################################################
##                                                                            ##
##                               fit_OU_trend   (v1.1)                        ##
##                                                                            ##
##  Created by Julien Clavel - 22-05-2019                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr).                  ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################


fit_t_env_ou <- function(phylo, data, env_data, error=NULL, model, method="Nelder-Mead", control=list(maxit = 20000),...){
  
  # options
  args = list(...)
  if(is.null(args[["upper"]])) upper <- Inf else upper <- args$upper
  if(is.null(args[["lower"]])) lower <- -Inf else lower <- args$lower
  if(is.null(args[["fixedRoot"]])) fixedRoot <- TRUE else fixedRoot <- args$fixedRoot
  if(is.null(args[["echo"]])){ echo <- TRUE }else{ echo <- args$echo}
  if(is.null(args[["param"]])){ param <- NULL }else{ param <- args$param}
  
  
  # parameters and objects needed for the computations
  nobs <- Ntip(phylo)
  C1 <- vcv.phylo(phylo)
  W <- matrix(0,ncol=1,nrow=nobs) # to tweak mvmorph function
  times <- diag(C1)
  root_age <-  max(node.depth.edgelength(phylo))
  
  # Max difference between tips and present day
  if(is.null(args[["maxdiff"]])){
    maxdiff<-0
  }else{
    maxdiff<-par$maxdiff
  }
  
  # reorder
  if(!is.null(names(data))) data <- data[phylo$tip.label] else warning("data and tips assumed to be aligned")
  
  if(is.null(env_data)){
    stop("Please provide a time-function or a time-serie dataset for the environmental changes; see ?fit_t_env_ou")
  }else if(!is.function(env_data)){
    
    # We check first that the climatic curve matches the phylogeny
    if(max(nodeHeights(phylo)) > max(env_data[,1])) stop("The environmental data does not cover the time span of the phylogeny: either enter data that covers the full time span or run analyses on younger clades")
    
    # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
    if(is.null(args[["df"]])){
      args$df <- smooth.spline(env_data[,1], env_data[,2])$df
    }
    spline_result <- sm.spline(env_data[,1],env_data[,2], df=args$df)
    env_func <- function(t){predict(spline_result,t)}
    
    # if we don't provide a time step in par we take the time steps of the dataset?
    t<-unique(env_data[,1])
    
    # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
    # the user can choose by specifying it in the "par" list
    
    if(is.null(args[["scale"]])){
      args$scale <- FALSE
    }
    
    # We build the interpolated smoothing spline function
    if(args$scale==FALSE){
      env_data<-splinefun(t,env_func(t))
    }else{
      curve_int<-env_func(t)
      curve_scaled=scale(curve_int,min(curve_int,na.rm=T),max(curve_int, na.rm=T)-min(curve_int,na.rm=T))
      env_data<-splinefun(t,curve_scaled)
    }
     
    # Otherwise we assume that the environmental function is given
  } 
  
  ## Parameterization
  if(missing(model)) model = NULL
  if(!is.function(model)){
    # use a linear function default model if no functions provided
    model<-function(x, env, param, theta0) theta0 + param[1]*env(x)
    if(is.null(args[["beta"]])) param <- c(0.1) else param <- args$beta
    if(echo==TRUE) message("Default model: theta(t) = theta_0 + beta * Env(t)")
  }

  # if not using default, then we should provide the starting values in the "param" argument
  if(is.null(param)) stop("You must provide starting values for your customized function")
  
  
  # Measurement error 
  nuisance <- FALSE
  if(!is.null(error)){
     # reorder the trait vector according to the tree
    if(!any(is.na(error))){
      if(is.null(names(error))){
        stop("You should provide a named vector for \"error\" ")
      }else{
        error<-error[phylo$tip.label]
      }
      nuisance <- TRUE # previous versions were error_param <- FALSE; to harmonize with fit_t_standard and fit_t_comp
    }else{
      # if NA is provided to "error", then we can estimate nuisance even if we don't have known measurement errors
      error <- rep(0, Ntip(phylo))
      nuisance <- TRUE
    }
  }
  
  # Define a function to find the OU expectation
  Expectation_OU <- function(time, theta_0, alpha, fun_exp, par, env){
    
    # Integral
    fun_ou_expectation <- function(x, stop_time, par, theta_0) { alpha*exp(alpha*(x-stop_time))*fun_exp((root_age+maxdiff)-x, env, par, theta_0) }
   
    # Numerical integration
    expectation_vector <- sapply(time, function(t){
       int <- try(integrate(f=fun_ou_expectation, lower=0, upper=t, stop_time=t, par=par, theta_0=theta_0, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent=TRUE)
      if(inherits(int ,'try-error')){
        warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
        integ <- NA_real_
      } else {
        integ <- int$value
      }
     #theta_0*exp(-alpha*t) + integ
      # TO CHECK but I think it's better to use as the start of the process X(0) = X(0)+beta*T(0) to take into account possible offset if the climatic curve at t=0 is not equal to 0
      exp(-alpha*t)*fun_exp((root_age+maxdiff), env, par=par, theta_0) + integ
    })
    return(expectation_vector)
  }
  
  # log-likelihood function
  llik <- function(par, phy, data, fun_ou, env){
    
    # parameters
    sigma2 = exp(par[1])
    alpha = exp(par[2])
    theta_0 = par[3]
    if(nuisance) extra_par = par[5:length(par)] else extra_par = par[4:length(par)]
    
    # Variance-Covariance
    if(fixedRoot){
       
        V<-.Call(C_panda_covar_ou_fixed, A=C1, alpha=alpha, sigma=sigma2)
     
        # add ME and nuisance variance
        if(!is.null(error) & nuisance==TRUE){
          nuisanceME = exp(par[4])
          diag(V) <- diag(V) + error^2 + nuisanceME
        }
        
    }else{
        
        V<-.Call(C_panda_covar_ou_random, A=C1, alpha=alpha, sigma=sigma2)
        
        # add ME and nuisance variance
        if(!is.null(error) & nuisance==TRUE){
          nuisanceME = exp(par[4])
          diag(V) <- diag(V) + error^2 + nuisanceME
          }
        
    }
    
    # Expectation
    E <- Expectation_OU(times, theta_0=theta_0, alpha=alpha, fun_exp=fun_ou, par=extra_par, env=env)
    # check if there were an error in the integration
    if(any(is.na(E))) return(list(logl=-1e7))
    
    # ll computation
    residuals <- (data - E)
    llik <- try(mvLL(V, residuals, method="rpf", param=list(D=W, estim=FALSE, mu=1)), silent = TRUE)
    
    # Try catch if there is an error
    if(inherits(llik ,'try-error')){
        warning("An error occured during the likelihood estimation. An infinite value was returned to continue the search")
        llik <- list()
        llik$logl <- 1e6
    }
    
    return(llik)
  }
  
  # startValues
  #sig_est <- mvLL(phylo, data, param=list(estim=TRUE))$sigma 
  if(nuisance){
    if(echo==TRUE) message("Finding best starting values...")
    #sig_start <- log(sig_est*c(0.2,0.5,0.8,0.9))
    #nuis_start <- log(sig_est*c(0.01, 0.05, 0.1))
    alpha_start <- log(log(2)/(max(node.depth.edgelength(phylo))/c(0.1,0.5,1.5,3,8)))
    sig_start <- log(Stationary_to_scatter(exp(alpha_start), var(data))) # condition the starting values on alpha?
    nuis_start <- log(exp(mean(sig_start))*c(0.01, 0.05, 0.1))
    theta_start <- mean(data)*c(0.5,0.8,1)
    start_val_ou <- list(sigma2=sig_start, alpha=alpha_start, theta_0=theta_start, mserr=nuis_start)
    
    mat_start <- expand.grid(c(start_val_ou ,param))
    start_val <- mat_start[which.min(apply(mat_start, 1, function(par) -llik(par, phy=phylo, data=data, fun_ou=model, env=env_data)$logl)),]#options
    names(start_val) = NULL
    
  }else{
    if(echo==TRUE) message("Finding best starting values...")
    # initialize the parameter search
    alpha_start <- log(log(2)/(max(node.depth.edgelength(phylo))/c(0.1,0.5,1.5,3,8)))
    sig_start <- log(Stationary_to_scatter(exp(alpha_start), var(data))) # condition the starting values on alpha?
    theta_start <- mean(data)*c(0.5,0.8,1)
    start_val_ou <- list(sigma2=sig_start, alpha=alpha_start, theta_0=theta_start)
    
    mat_start <- expand.grid(c(start_val_ou ,param))
    start_val <- mat_start[which.min(apply(mat_start, 1, function(par) -llik(par, phy=phylo, data=data, fun_ou=model, env=env_data)$logl)),]#options
    names(start_val) = NULL
  } 
  
  # optimization
  if(echo==TRUE) message("Start optimization (and numerical integration). Please wait...")
  # if(method=="spg"){
  #     require(BB)
  # estimModel <- spg(unlist(start_val),
  #     fn = function(par) -llik(par, phy=phylo, data=data, fun_ou=env_data)$logl,
  #     method=3,
  #     upper=upper,
  #     lower=lower,
  #     control=control)
  # }else{
  estimModel <- optim(start_val,
                      fn = function(par) -llik(par, phy=phylo, data=data, fun_ou=model, env=env_data)$logl,
                      method=method,
                      upper=upper,
                      lower=lower,
                      hessian=TRUE,
                      control=control)
  #}
  
  # Done retrieve param
  if(echo==TRUE) message("Done. retrieve parameters and results...")
  param = estimModel$par
  LL = -estimModel$value
  nparam = length(param)
  
  # parameters
  sigma2 = exp(param[1])
  alpha = exp(param[2])
  theta_0 = param[3]
  if(nuisance) estimated_ME = exp(param[4]) else estimated_ME = NULL
  if(nuisance) custom = param[5:length(param)] else custom = param[4:length(param)]
  
  if(nuisance){
  list_par = list(sigma2=sigma2,
                        alpha=alpha,
                        theta_0=theta_0,
                        par=custom,
                        nuisance = estimated_ME,
                        root_implied=model(0, env_data, custom, theta_0))
  }else{
      list_par = list(sigma2=sigma2,
                            alpha=alpha,
                            theta_0=theta_0,
                            par=custom,
                            root_implied=model(0, env_data, custom, theta_0))
  }
  
  # AIC
  AIC = -2*LL+2*nparam
  # AIC corrected
  AICc = AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) 
  
  # Check the curvature of the likelihood surface from the Hessian to assess the convergence (does not guarantee that it's a local optima...)
  # if(method=="spg"){
  #     library(numDeriv)
  #     hess<-eigen(hessian(function(par) -llik(par, phy=phylo, data=data, fun_ou=env_data)$logl, param))$values
  # }else{
      hess<-eigen(estimModel$hessian)$values
  # }
  if(any(hess<0)) hess.value<-1 else hess.value<-0
  
  # Return the results
  results <- list(LH=LL, aic=AIC, aicc=AICc, free.parameters=nparam, param=list_par, root=list_par$root_implied, opt=estimModel, convergence=estimModel$convergence, hess.value=hess.value, env_func=env_data, model=model,
                  nuisance=estimated_ME, tot_time=root_age+maxdiff)

  class(results) <- "fit_t.env.ou"
  return(results)
  #invisible(results)
}


# Print the results
print.fit_t.env.ou<-function(x,...){
    cat("OU model with time-dependent optimum","\n")
    cat("AIC :",x$aic,"\n")
    cat("AICc:",x$aicc,"\n")
    cat("Log-Likelihood:",x$LH,"\n")
    cat("______________________","\n")
    cat("Parameters:","\n")
    print(unlist(x$param), digits=3)
    cat("______________________","\n")
    if(x$convergence==0){cat("Succesful convergence","\n")}else{cat("Convergence has not been reached","\n")}
    if(x$hess.value==0 & x$convergence==0){cat("A reliable solution has been reached","\n")}else{cat("Unreliable solution (Likelihood at a saddle point)","\n")}
}

Stationary_to_scatter <- function(alpha, stationary){
  return(alpha*stationary + stationary*(alpha))
}

# # Test
# n=100
# tree <- rtree(n)
# tree$edge.length <- tree$edge.length * 60/ max(nodeHeights(tree))
# 
# # define a function for the optimum theta(s)
# fun_temp <- function(x, par)  par[1]*d18(max(nodeHeights(tree))-x) # here it's max-x because the time start at present
# 
# # simulate random data
# dat <- rTraitCont(tree)
# 
# # fit the model
# fit_OU_trend(tree, dat, fun=fun_temp, param=c(1))
