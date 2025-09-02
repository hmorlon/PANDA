################################################################################
##                                                                            ##
##                               fit_threshold  (v1.0)                        ##
##                                                                            ##
##  Created by Julien Clavel - 01-09-2025                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##                                                                            ##
##                                                                            ##
################################################################################

## -------------------------------------------------------------------------------------- ##
#   Environmental Threshold Model for binary character - Julien Clavel 2023 - julien.clavel@univ-lyon1.fr/julien.clavel@cnrs.fr
# Uses Hiscott et al. 2016 (GBE) pruning algorithm with numerical integration to approximate the log-likelihood (adaptated from Mathlab code).
# The algorithm is used to compute ML estimate of Felsenstein' threshold model with the climatic model described in Clavel & Morlon 2017 (PNAS), the BM, EB, OU and lambda models
# "phylo" is an object of class phylo (see ?ape)
# data is vector of a discrete (binary) trait (order should match the order in the tree)
# N = is the number of knots used to approximate the integral (through numerical integration) in the calculation of the partial likelihoods. The higher, the better, at the expense of computation time.
# env_data is an environmental curve object. See ?fit_t_env in RPANDA package for further details
## -------------------------------------------------------------------------------------- ##

# NOTE: code is quite slow, should implement the numerical integration in C or find some alternatives

# Workhorse function for the various Threshold models
fit_t_thresh <- function(phylo, data, model=c("BM","OU","EB","lambda","Clim"), N=200, env_data=NULL, ...){
    
    # Reorder the tree => tree to phylo !!!
    if(!inherits(phylo,"postorder")) phylo<-reorder.phylo(phylo,"postorder")
    
    # reorder the trait vector according to the tree
    if(is.null(names(data))){
        stop("You should provide a named vector for \"data\" ")
    }else{
        data<-data[phylo$tip.label]
    }
    
    # TODO: check data type, should be a two-state character
    
    # check model
    model = match.arg(model[1],c("BM","OU","EB","lambda","Clim"))
    
    ## Use ellipsis for param arguments
    para<-list(...)
    tot_time = 0; # only used for the Clim model
    
    # switch between methods
    switch(model,
    "BM"={
        fun_to_optim <- function(par){
          ll <- thresholdLikelihoodGauss(phylo, data, N, par[1])
          return(-ll)
        }
        
        # parameter to optimize: only the threshold value
        start_guess = 0
        lower = -Inf
        upper = Inf
    },
    "EB"={
        fun_to_optim <- function(par){
          if(abs(par[1])>.Machine$double.eps) phy2 <- transform_EB(phylo, beta=par[1], sigmasq=1) else phy2<-phylo
          ll <- thresholdLikelihoodGauss(phy2, data, N, par[2])
          if(!is.finite(ll)) return(1e6)
          return(-ll)
        }
        
        maxHeight <- max(node.depth.edgelength(phylo))
        start_guess <-  c(-log(2)/(maxHeight/1.5), 0)
        lower <- c(log(10^-5)/maxHeight, -Inf)
        upper <- c(0, Inf)
    },
    "OU"={
        if(!is.ultrametric(phylo)) stop("This model is currently not suported for non-ultrametric trees")
        if(is.null(para[["alpha"]])) start_guess =  log(2)/(max(node.depth.edgelength(phylo))/1.5) else start_guess = para$alpha
        
        fun_to_optim <- function(par){
          phy2 <- .phyOU(phylo, par[1])
          ll <- thresholdLikelihoodGauss(phy2, data, N, par[2])
          if(!is.finite(ll)) return(1e6)
          return(-ll)
        }
        
        # define default upper bound for alpha parameter search (~25 phylogenetic half-lives)
        start_guess <- c(start_guess, 0)
        upper <- c(log(2)/(max(node.depth.edgelength(phylo))/25), Inf)
        lower <- c(1e-10,-Inf)
        
    },
    "lambda"={
        
        fun_to_optim <- function(par){
          phy2 <- lambda_transform(phylo, par=par[1])
          ll <- thresholdLikelihoodGauss(phy2, data, N, par[2])
          return(-ll)
        }
        
        # define the bounds for the parameter search
        upper <- c(1, Inf)
        lower <- c(1e-6, -Inf)
        
        # start_guess
        start_guess = c(0.5, 0)
    },
    "Clim"={
        # Control for the integrate function: number of subdivisions
        if(is.null(para[["subdivisions"]])){
          subdivisions<-200L
        }else{
          subdivisions<-para$subdivisions
        }
        
        # Max difference between tips and present day
        if(is.null(para[["maxdiff"]])){
          maxdiff<-0
        }else{
          maxdiff<-para$maxdiff
        }
        
        # Compute the branching times
        if(is.ultrametric(phylo)){
          times<-branching.times(phylo)
        }else{
          # Use "phytools" called by mvMORPH
          
          times<-max(nodeHeights(phylo))-nodeHeights(phylo)[match(1:phylo$Nnode+Ntip(phylo),phylo$edge[,1]),1]
          names(times)<-1:(phylo$Nnode+Ntip(phylo))
        }
        
        # Root age
        tot_time<-max(times)
        
        # Set the root to zero
        times<-tot_time-times
        
        # Check if the climatic function is provided
        if(is.null(env_data)){
          stop("Please provide a time-function or a time-serie dataset for the environmental changes; see ?fit_t_env")
        }else if(!is.function(env_data)){
          
          # We check first that the climatic curve matches the phylogeny
          if(max(nodeHeights(phylo)) > max(env_data[,1])) stop("The environmental data does not cover the time span of the phylogeny: either enter data that covers the full time span or run analyses on younger clades")
          
          # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
          if(is.null(para[["df"]])){
            para$df <- smooth.spline(env_data[,1], env_data[,2])$df
          }
          spline_result <- sm.spline(env_data[,1],env_data[,2], df=para$df)
          env_func <- function(t){predict(spline_result,t)}
          
          # if we don't provide a time step in par we take the time steps of the dataset?
          t<-unique(env_data[,1])
          
          # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
          # the user can choose by specifying it in the "par" list
          
          if(is.null(para[["scale"]])){
            para$scale <- FALSE
          }
          
          # We build the interpolated smoothing spline function
          if(para$scale==FALSE){
            env_data<-splinefun(t,env_func(t))
          }else{
            curve_int<-env_func(t)
            curve_scaled=scale(curve_int,min(curve_int,na.rm=T),max(curve_int, na.rm=T)-min(curve_int,na.rm=T))
            env_data<-splinefun(t,curve_scaled)
          }
          
          # Otherwise we assume that the environnemental function is given
        }
        
        # wrapper to the main Clim-transform function
        transClim <- function(tree, beta, sigma=1){
            transformClim_thresh(tree, param=c(sigma,beta), fun=env_data, times=times, check=FALSE, mtot=tot_time, subdivisions=subdivisions, maxdiff=maxdiff, model="EnvExp")
        }
        
        # Function to optimize
        fun_to_optim <- function(par){
          phy2 <- transClim(phylo, beta=par[1])
          if(any(!is.finite(phy2$edge.length))) return(1e6)
          ll <- thresholdLikelihoodGauss(phy2, data, N, par[2])
          if(!is.finite(ll)) return(1e6)
          return(-ll)
        }
        
        # define bounds for parameters search
        upper = c(Inf,Inf)
        lower = c(-Inf,-Inf)
        start_guess = c(0., 0)
        
    })
    
    
    # Optimization of the threshold model log-likelihood
    opt <- optim(par=start_guess, fn = fun_to_optim,  method="L-BFGS-B", lower=lower, upper=upper)
    
    # results
    nparam = length(opt$par)
    nobs = length(data)
    AIC = 2*opt$value+2*nparam
    AICc = AIC+((2*nparam*(nparam+1))/(nobs-nparam-1))
    
    # estimated parameters
    switch(model,
    "BM"={parameters <- data.frame(sigma=1, mu=opt$par[1])},
    "lambda"={parameters <- data.frame(lambda=opt$par[1], mu=opt$par[2])},
    "OU"={parameters <- data.frame(alpha=opt$par[1], mu=opt$par[2])},
    "EB"={parameters <- data.frame(a=opt$par[1], mu=opt$par[2])},
    "Clim"={parameters <- data.frame(beta=opt$par[1], mu=opt$par[2])})
    
    # return the results
    results <- list(LH=-opt$value,
                      free.parameters=nparam,
                      param=parameters,
                      aic=AIC, aicc=AICc, opt=opt, model=model, convergence=opt$convergence, env_func=env_data, tot_time=tot_time)
    class(results) = "threshML"
    return(results)
}



# pruning algorithm with numerical integration
thresholdLikelihoodGauss <- function(phylo, data, N=200, mu=0){
  
  if(!inherits(phylo,"postorder")) phylo = reorder.phylo(phylo, "postorder")
  
  # parameters
  n = Ntip(phylo)
  nnodes = 2*n-1
  nchar = length(unique(data))
  C = model.matrix(~as.factor(data)+0)
  
  # vecteurs, matrices
  L = numeric(nnodes)
  U = numeric(nnodes)
  X = matrix(0, nrow=nnodes, ncol=N+1)
  x = (0:N)/N # for linear interpolation
  
  # truncation error for each individual integral
  alpha = n*N^(-4)/(2*(n-1))
  
  # populate X value for the root? given or fixed
  X[n+1,] = mu
  
  # populate L, U and X
  for(i in (n+2):nnodes){
    L[i] = qnorm(alpha, mu, sqrt(sum(phylo$edge.length[hit_tip(i,phylo)])))
    U[i] = qnorm(1-alpha, mu, sqrt(sum(phylo$edge.length[hit_tip(i,phylo)])))
    X[i,] = (U[i] - L[i])*x+L[i]
  }
  
  # parameters
  logL = 0 # log of the product of the likelihoods
  like = numeric(nchar)
  pair = seq(from=1, to=(2*n-4), by=2) # to (2*N-2) but (2*N-4) avoid updating the ancestral state
  edge1 = phylo$edge[,1]
  edge2 = phylo$edge[,2]
  
  for(ch in 1:nchar){
    z = C[,ch] # traits values observed at the leaf for characters in C
    F = matrix(1, nnodes, N+1)  #F(i,x) = likelihood of the data at or directly
    #below node i, given a liability of x at i's
    #parent node.
    
    # tips values
    for(i in 1:n){
      tip_i = which(edge2==i)
      j = edge1[tip_i] # parent of i
      
      if(z[i]==1){
        # integrates over positive liabilities
        F[i,] = 1 - pnorm(0, X[j,], sqrt(phylo$edge.length[tip_i]))
      }else if(z[i]==0){
        # integrates over negative liabilities
        F[i,] = pnorm(0, X[j,], sqrt(phylo$edge.length[tip_i]))
      }
    }
    
    # pruning down the tree (all internal nodes excluding the root)
    for(h in pair){
      i = edge1[h]              # node "i"
      i_index = which(edge2==i) # index of node "i"
      j = edge1[i_index]        # parent of node "i"
      u = edge2[h]              # left descendent of node "i"
      v = edge2[h+1]            # right descendent of node "i"
      
      if(phylo$edge.length[i_index]<=.Machine$double.eps){
        F[i,] = F[u,]*F[v,]
      }else{
        # calculate the integral using a kernel (e.g., simpson, gaussian)
        weights = buildWeightsSet(X[i,],X[j,], phylo$edge.length[i_index])
        F[i,] <- tcrossprod(F[u,]*F[v,], weights) # that is: weights%*%(F[u,]*F[v,]); must check if it's faster than using 'for' loops
      }
    }
    
    # log-likelihood at the root (assumes it's a fixed value)
    u = which(edge1==n+1)
    # We have L_r = U_r, so X_r[k] is constant in k.
    like[ch] = F[edge2[u[1]],1]*F[edge2[u[2]],1]
    
    # check first that the value is safe?
    if(!is.finite(like[ch])){
      warning("NaN were produced - there were an issue in the likelihood computation")
      return(-Inf)
    }
    
    if(like[ch]>0){
      logL = logL + log(like[ch])
    }else{
      logL = -Inf
    }
  }
  
  return(logL)
}


# recover the path from root to tips
hit_tip <- function(tip, tree){
  
  # parameters
  N <- Ntip(tree)
  edge1 = tree$edge[,1]
  edge2 = tree$edge[,2]
  i=2
  list_ind <- list()
  list_ind[[1]] <- tip
  while(tip!=(N+1)){
    list_ind[i] <- tip <- edge1[tree$edge[,2]==tip]
    i = i+1
  }
  
  node_path<-unlist(list_ind)
  index = which(tree$edge[,2]%in%node_path[-length(node_path)])
  
  return(index)
}



# compute weights for the Gaussian kernel quadrature
buildWeightsSet <- function(X, mu, sig2){
  #weights calculated for approximating \int phi(x,mu,sig2) f(x) dx
  N = length(X)-1;
  w = matrix(0,length(X),length(X)); #weight array
  h = X[2] - X[1]; # spacing between two consecutive points.
  f = matrix(0, nrow=3, ncol=(N/2)+1); #matrix of antiderivatives used to calculate the weights.
  sig = sqrt(sig2);
  index = seq(from=1,to=(N+1),by=2)
  x = matrix(X[index], nrow=N+1, ncol=N/2+1, byrow = TRUE)
  MU = matrix(t(mu), ncol=N/2+1, nrow = length(X))
  #f(1,i) = antiderivative of phi(x,mu,sig2) at x = X(2*i - 1)
  f1 = .5*erf((x-MU)/(sqrt(2*sig2)));
  #f(2,i) = antiderivative of x*phi(x,mu,sig2) at x = X(2*i - 1)
  f2 = -sig2*dnorm(x, MU, sig) + MU*f1;
  #f(3,i) = antiderivative of x^2*phi(x,mu,sig2) at x = X(2*i - 1)
  f3 = (-2*sig*((x+MU)*exp(-((x-MU)^2)/(2*sig2))) + 2*sqrt(2*pi)*(MU^2 + sig2)*f1)/(2*sqrt(2*pi));
  x = x[,-ncol(x)];
  x2 = x^2;
  F1 = f1[,-1] - f1[,-ncol(f1)]# Note that \int f(x) dx [a, b] = F(b) - F(a) where F is the antiderivative of f
  F2 = f2[,-1] - f2[,-ncol(f2)]
  F3 = f3[,-1] - f3[,-ncol(f3)]
  z = matrix(0, ncol=1, nrow=N+1)
  w[,seq(from=1, to=ncol(w), by=2)] = (.5*(((2*h^2)+3*h*cbind(x,0)+cbind(x2,0))*cbind(F1,0) + (cbind(0,x)*(cbind(0,x+h))*cbind(0,F1))) - .5*(((3*h)+(2*cbind(x,0)))*cbind(F2,0) + ((2*cbind(0,x))+h)*cbind(0,F2)) + .5*(cbind(F3,0) + cbind(0,F3)))/(h^2);
  w[,seq(from=2, to=ncol(w), by=2)] = (-((x2+2*h*x)*F1)+(2*(x+h)*F2)-F3)/(h^2);
  return(w)
}

# the so-called 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# Function for the Pagel's lambda model. Stretches the branch lengths
lambda_transform <- function(phylo, par){
    
    n <- Ntip(phylo)
    parent <- phylo$edge[,1]
    descendent <- phylo$edge[,2]
    extern <- (descendent <= n)
    
    # Pagel's lambda tree transformation
    if(par!=1) {
        root2tipDist <- node.depth.edgelength(phylo)[1:n] # for non-ultrametric trees. The 'up' limit should be exactly 1 to avoid singularity issues
        phylo$edge.length <- phylo$edge.length * par
        phylo$edge.length[extern] <- phylo$edge.length[extern] + (root2tipDist * (1-par))
    }
    
    return(phylo)
}

# Wrapper to the .CLIMtransform function used in fit_t_env
transformClim_thresh<-function(phylo, model=c("EnvExp", "EnvLin"), ...){
  
  # Use mvMORPH for computing the log-likelihood
  #require(mvMORPH)
  
  ## Parameterization
  par<-list(...)
  
  # Default model
  if(!is.function(model)){
    model<-model[1]
  }
  
  # Number of tips
  tips <- length(phylo$tip.label)
  
  # Check if the climatic function is provided
  if(is.null(par[["fun"]])){
    stop("Please provide a time-function")
  }else if(!is.function(par$fun)){
    
    # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
    if(is.null(par[["df"]])){
      par$df <- smooth.spline(par$fun[,1], par$fun[,2])$df
    }
    spline_result <- sm.spline(par$fun[,1],par$fun[,2], df=par$df)
    env_func <- function(t){predict(spline_result,t)}
    
    # if we don't provide a time step in par we take the time steps of the dataset?
    t<-unique(par$fun[,1])
    
    # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
    # the user can choose by specifying it in the "par" list
    
    if(is.null(par[["scale"]])){
      # We build the interpolated smoothing spline function
      par$fun<-splinefun(t,env_func(t))
    }else{
      curve_int<-env_func(t)
      curve_scaled=scale(curve_int,min(curve_int,na.rm=T),max(curve_int, na.rm=T)-min(curve_int,na.rm=T))
      par$fun<-splinefun(t,curve_scaled)
    }
    
    # Otherwise we assume that the environnemental function is given
  }
  
  # Check if the branching time is provided
  if(is.null(par[["times"]])){
    warning("The branching time for the \"phylo\" object was not provided by the user")
    par$times<-branching.times(phylo)
    # root age
    mtot<-max(par$times)
    # Set the root to zero
    par$times<-max(par$times)-par$times
  }else{
    # root age
    mtot<-par$mtot
  }
  
  # Check if the difference between tip and present day.
  if(is.null(par[["maxdiff"]])){
    maxdiff<-0
  }else{
    maxdiff<-par$maxdiff
  }
  
  
  # Check if the tree is in prunning-wise order
  if(is.null(par[["check"]])){
    par$check<-TRUE
  }
  
  # Check if there is polytomies
  if (phylo$Nnode != tips - 1) {
    stop("You can't use this function with polytomies, transform the tree using \"multi2di\" function first")
  }
  
  
  # Control for the integrate function: number of subdivisions
  if(is.null(par[["subdivisions"]])){
    subdivisions<-200L
  }else{
    subdivisions<-par$subdivisions
  }
  
  ## Transform the tree and return the log-likelihood
  
  # Check the parameters
  if(is.null(par[["param"]]))  {
    stop("Please provide parameters values for \"sig2\" and \"beta\" ")
  }
  
  # Sigma is not provided but analytically computed instead
  phylo <- .CLIMtransform(phylo, param=par$param, mtot=mtot, times=par$times, funEnv=par$fun, model=model, tips=tips, subdivisions=subdivisions, maxdiff=maxdiff)
  
  return(phylo)
}


# Print the results
print.threshML<-function(x,...){
    cat("Threshold model",x$model,"\n")
    cat("AIC :",x$aic,"\n")
    cat("AICc:",x$aicc,"\n")
    cat("Log-Likelihood:",x$LH,"\n")
    cat("______________________","\n")
    cat("Parameters:","\n")
    print(unlist(x$param), digits=3)
    cat("______________________","\n")
    if(x$opt$convergence==0){cat("Succesful convergence","\n")}else{cat("Convergence has not been reached","\n")}
    #if(x$hess.value==0 & x$opt$convergence==0){cat("A reliable solution has been reached","\n")}else{cat("Unreliable solution (Likelihood at a saddle point)","\n")}
}
