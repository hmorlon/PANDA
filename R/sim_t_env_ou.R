##--------------------------------------------------------------------------------##
##                                                                                ##
##                    Simul trait OU Climatic (recursive)                         ##
##                                                                                ##
##--------------------------------------------------------------------------------##

sim_t_env_ou<-function(phylo, param, env_data, model, step=0.01, plot=FALSE, sigma, alpha, theta0, ...){
  
  # options
  arg <- list(...)
  if(missing(phylo)) stop("You must provide a phylogenetic tree!!")
  if(missing(param)) stop("You must provide a vector of parameters or a model fit object in the \"param\" argument")
  
  # model fit and options
  if(inherits(param,"fit_t.env.ou")){
    fun <- param$model
    env_data <- param$env_func
    root.value <- param$root
    theta0 <- param$param$theta_0
    sigma <- sqrt(param$param$sigma2)
    alpha <- param$param$alpha
    param <- param$param$par # will erase the previous param argument !!!
  }else{
    if(missing(env_data)) stop("You must provide environmental data!!")
    if(missing(sigma)) stop("You must provide sigma parameter!!")
    if(missing(alpha)) stop("You must provide alpha parameter!!")
    if(missing(theta0)) stop("You must provide theta0 parameter!!")
    if(missing(model)) model = NULL
    
    # the default function model
    if(!is.function(model)){
      fun <- function(x, env_data, param, theta0){theta0 + param[1]*env_data(x)}
      message("Default model: theta(t) = theta_0 + beta * Env(t)")
    }else{
      fun <- model
    }
  }
  
  # Check if the climatic function is provided
  if(!is.function(env_data)){
    
    # env_data is a dataframe with two columns: 1) is time; 2) is datapoints
    if(is.null(arg[["df"]])){
      arg$df <- smooth.spline(env_data[,1], env_data[,2])$df
    }
    spline_result <- sm.spline(env_data[,1],env_data[,2], df=arg$df)
    env_func <- function(t){predict(spline_result,t)}
    
    # if we don't provide a time step in par we take the time steps of the dataset?
    t<-unique(env_data[,1])
    
    # If the ClimLin model is used we should standardize the environmental curve: i.e. sigma is the maximum rate value.
    # the user can choose by specifying it in the "par" list
    
    if(is.null(arg[["scale"]])){
      # We build the interpolated smoothing spline function
      env_data<-splinefun(t,env_func(t))
    }else{
      curve_int<-env_func(t)
      curve_scaled=scale(curve_int,min(curve_int,na.rm=TRUE),max(curve_int, na.rm=TRUE)-min(curve_int,na.rm=TRUE))
      env_data<-splinefun(t,curve_scaled)
    }
    
    # Otherwise we assume that the environnemental function is given
  }
  
  ## First define the function for simulating the SDEs
  
  ou_model<- function(startx, nwalks=50, theta, sigma, alpha, dt_time, dt, remind=FALSE, dt2=0, param, env, theta0)
  {
    x <- numeric(nwalks)
    x[1] <- startx
    for (i in 1:(nwalks-1))
    {
      if(remind==TRUE && i==(nwalks-1)){
        x[i+1] <- x[i] + dt2*alpha*(theta(dt_time[i], env, param, theta0)-x[i]) + rnorm(1, sd=sigma*sqrt(dt2))
      }else{
        x[i+1] <- x[i] + dt*alpha*(theta(dt_time[i], env, param, theta0)-x[i]) + rnorm(1, sd=sigma*sqrt(dt))
      }
      
    }
    return(x)
  }
  
  ## Run the process
  SimT<-function(start, n, dt, sig=1, alpha, age, fun, rootage, param, env, theta0){
    if(n%%dt==0){
      integ<-n%/%dt
      if(integ!=0){
        ageval<-rootage-(cumsum(rep(dt,integ))+age)
        
        # run OU process on the branch
        process <- ou_model(start, nwalks=(integ+1), theta=fun, sigma=sig, alpha=alpha, dt_time=ageval,  dt=dt, param=param, env=env, theta0=theta0)
        
        
      }else{
        ageval<-rootage-((n%%dt)+age)
        process <- c(start,(n%%dt)*alpha*(fun(ageval, env, param, theta0)-start) + rnorm(1, sd=sig*sqrt((n%%dt))) )
      }
      
    }else{
      integ<-(n%/%dt)
      if(integ!=0){
        ageval<-rootage-(cumsum(c(rep(dt,integ),n%%dt))+age)
        process <- ou_model(start, nwalks=length(ageval), theta=fun, sigma=sig, alpha=alpha, dt_time=ageval, dt=dt, remind=TRUE, dt2=n%%dt, param=param, env=env, theta0=theta0) # should include also the integration of n%%dt
      }else{
        ageval<-rootage-((n%%dt)+age)
        process <- c(start,(n%%dt)*alpha*(fun(ageval, env, param, theta0)-start) + rnorm(1, sd=sig*sqrt(n%%dt)))
      }
    }
    return(process)
  }
  
  ## discretize time
  discT<-function(n,dt,age,rootage){
    if(n%%dt==0){
      integ<-n%/%dt
      if(integ!=0){
        sigmaval<-rootage-(cumsum( c(0,rep(dt,integ)) )+age)
      }else{
        sigmaval<-rootage-(c(0,n%%dt)+age)
      }
      
    }else{
      integ<-(n%/%dt)
      if(integ!=0){
        sigmaval<-rootage-(cumsum(c(rep(dt,integ),n%%dt))+age)#ne pas standardiser par sig puisque juste pour le plot...
      }else{
        sigmaval<-rootage-(c(0,n%%dt)+age)
      }
    }
    return(sigmaval)
  }
  
  
  tr<-reorder(phylo,"cladewise")
  Ages<-nodeHeights(tr)[,1]
  rootage<-max(nodeHeights(phylo))
  X<-T<-list()
  N<-length(tr$tip)
  for(i in 1:nrow(tr$edge)){
    if(tr$edge[i,1]==(N+1)){
      # temp root age or just theta?
      # if(theta0){
      #   X[[i]]<-SimT(start=anc, n=tr$edge.length[i], dt=step, sig=sigma, alpha=alpha, age=Ages[i], env=env_data, fun=fun, rootage=rootage, param=param, theta0=theta0) # N+1 = element ?? la racine; anc =state at the root
      # }else{
      #   ## Assume here that env_data(0) is present days, so we need to specify rootage instead
        X[[i]]<-SimT(start=fun(rootage, env_data, param, theta0), n=tr$edge.length[i], dt=step, sig=sigma, alpha=alpha, age=Ages[i], env=env_data, rootage=rootage, fun=fun, param=param, theta0=theta0)
      # }# N+1 = element ?? la racine;
    }else{
      parent<-match(tr$edge[i,1],tr$edge[,2]) # recherche le noeud parent de l'element i (peut par ex. etre N+1 si deuxieme noeud sur l'arbre). renvoie la ligne dans edge
      X[[i]]<-SimT(start=X[[parent]][length(X[[parent]])],n=tr$edge.length[i],dt=step, sig=sigma, alpha=alpha, age=Ages[i], env=env_data, rootage=rootage, fun=fun, param=param, theta0=theta0)
    }
    T[[i]] <- discT(tr$edge.length[i], dt=step, age=Ages[i], rootage)
  }
  ## Plot function
  bm<-function(X,T,...){
    minX<-min(sapply(X,min))
    maxX<-max(sapply(X,max))
    plot(T[[1]][1],X[[1]][1],ylim=c(minX,maxX),xlim=c(0,max(sapply(T,max))),ylab="phenotype",xlab="time")
    for(i in 1:length(X)) lines(T[[i]],X[[i]])
    if(hasArg(colors)) cols<-list(...)$colors
    else cols<-"black"
    return(cols)
  }
  if(plot==TRUE)  cols<-bm(X,T)
  # retourne les valeurs
  ind<-sapply(1:N,function(x){which(tr$edge[,2]==x)})
  val<-sapply(1:N,function(x){X[[ind[x]]][length(X[[ind[x]]])] } )
  names(val) = tr$tip.label
  return(val)
}
