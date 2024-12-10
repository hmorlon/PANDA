.Psi_in_past <- function(s,t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=0)
{
  if ((cst.lamb==TRUE) & (cst.mu==TRUE)){
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- exp(r*(t-s))*(abs(1+(lamb*(exp(r*t)-exp(r*s)))/(r/f+lamb*(exp(r*s)-1))))^(-2)
    return(res)
  }
  
  ####### exponential dependencies ########
  
  if ((cst.lamb==TRUE) & (expo.mu==TRUE)){
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
    return(res)
  }
  
  if ((expo.lamb==TRUE) & (cst.mu==TRUE)){
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
    return(res)
  }
  
  if ((expo.lamb==TRUE) & (expo.mu==TRUE)){
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0/beta*(exp(beta*y)-exp(beta*x))}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
    return(res)
  }
  
  ####### other dependencies ########
  
  if (dt==0){
    # Compute using R integration functions
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int <- function(x,y){.Integrate(Vectorize(r),x,y,stop.on.error=FALSE)}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)}
    rst <- r.int(s,t)
    rist <- r.int.int(s,t)
    ri0s <- r.int.int(0,s)
    res <- exp(rst)*(abs(1+rist/(1/f+ri0s)))^(-2)
    return(res)
    
  }else{
    
    Nintervals <- 1 + as.integer(t/dt)
    X <- seq(0, t, length.out = Nintervals + 1)
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int <- cumsum(sapply(X, r)) * t / Nintervals
    r.int.0 <- function(y){exp(r.int[1 + as.integer( y * Nintervals / t)]) * f.lamb(y)}
    r.int.int.tab <- cumsum(sapply(X, r.int.0)) * t / Nintervals
    r.int.int <- function(x,y){
      indy <- 1 + as.integer( y * Nintervals / t)
      indx <- 1 + as.integer( x * Nintervals / t)
      value <- r.int.int.tab[indy] - r.int.int.tab[indx]
      return(value)
    }
    rst <- r.int[1 + Nintervals] - r.int[1 + as.integer( s * Nintervals / t)]
    rist <- r.int.int(s,t)
    ri0s <- r.int.int(0,s)
    res <- exp(rst)*(abs(1+rist/(1/f+ri0s)))^(-2)
    return(res)
  }
}


.Phi_in_past <- function(t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=0)
{
  
  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- 1-r*exp(r*t)/(r/f+lamb*(exp(r*t)-1))
    return(res)
  }
  
  if ((cst.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }
  
  if ((expo.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }
  
  if ((expo.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0/beta*(exp(beta*y)-exp(beta*x))}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }
  
  else
  {
    if (dt==0)
    {
      r <- function(t){f.lamb(t)-f.mu(t)}
      r.int <- function(x,y){.Integrate(Vectorize(r),x,y,stop.on.error=FALSE)}
      r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
      r.int.int <- function(x,y){.Integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)}
      rit <- r.int(0,t)
      ri0t <- r.int.int(0,t)
      res <- 1.0 - exp(rit)/(1/f+ri0t)
      return(res)
    }
    else
    {
      s <- 0
      t <- t
      Nintervals <- 1 + as.integer((t-s)/dt)
      X <- seq(s, t, length.out = Nintervals + 1)
      r <- function(t){f.lamb(t)-f.mu(t)}
      r.int <- cumsum(sapply(X, r)) * (t - s) / Nintervals
      r.int.0 <- function(y){exp(r.int[1 + as.integer( (y - s) * Nintervals / (t - s))]) * f.lamb(y)}
      r.int.int.tab <- cumsum(sapply(X, r.int.0)) * (t - s) / Nintervals
      r.int.int <- function(x,y){
        indy <- 1 + as.integer( (y - s) * Nintervals / (t - s))
        indx <- 1 + as.integer( (x - s) * Nintervals / (t - s))
        value <- r.int.int.tab[indy] - r.int.int.tab[indx]
        return(value)
      }
      rit <- r.int[1 + Nintervals]
      ri0t <- r.int.int(0,t)
      res <- 1.0 - exp(rit)/(1/f+ri0t)
      return(res)
    }
  }
}




likelihood_bd_in_past <- function(phylo,tot_time,time_stop,f.lamb,f.mu,desc, tot_desc,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=0,cond="crown"){
  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"")
  
  f <- desc/tot_desc
  
  nbtips <- Ntip(phylo)
  log_indLikelihood <- c()
  from_past <- cbind(phylo$edge,node.age(phylo)$ages)
  ages <- rbind(from_past[,2:3],c(nbtips+1,0))
  ages <- ages[order(ages[,1]),]
  age <- max(ages[,2])
  
  for (j in 1:(nbtips-1)){
    node <- (nbtips+j)
    edges <- phylo$edge[phylo$edge[,1]==node,]
    tj <- age-ages[edges[1,1],2]
    Psi_timevar_errap_tj <- .Psi_in_past(time_stop,tj,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)
    log_lik_tj <- log(f.lamb(tj)) + log(Psi_timevar_errap_tj)
    log_indLikelihood <- c(log_indLikelihood,log_lik_tj)
  }
  log_indLikelihood <- c(log_indLikelihood,log(.Psi_in_past(time_stop,tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)))
  
  log_data_lik <- sum(log_indLikelihood)+desc*log(f)
  
  log_data_lik <- log_data_lik + nbtips*log(1-.Phi_in_past(time_stop,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt))
  
  
  if (cond==FALSE){
    log_final_lik <- log_data_lik
  }
  
  else if (cond=="stem"){
    Phi <- .Phi_in_past(tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)
    log_final_lik <- log_data_lik - log(1-Phi)
  }
  
  else if (cond=="crown"){
    Phi <- .Phi_in_past(tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt)
    log_final_lik <- log_data_lik - log(f.lamb(tot_time)) - 2*log(1-Phi)
  }
  
  return(log_final_lik)
}


fit_bd_in_past <- function (phylo, tot_time, time_stop, f.lamb, f.mu, lamb_par, mu_par, desc, tot_desc,
                            meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE,
                            expo.lamb=FALSE, expo.mu=FALSE, fix.mu=FALSE,
                            dt=0, cond="crown") {
  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"")
  
  
  nobs <- Ntip(phylo)
  
  phylo$edge.length[phylo$edge[,2] %in% 1:nobs] <- phylo$edge.length[phylo$edge[,2] %in% 1:nobs] + time_stop
  
  if (fix.mu==FALSE){
    
    init <- c(lamb_par,mu_par)
    p <- length(init)
    optimLH <- function(init){
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1+length(lamb_par)):length(init)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
      
      LH <- likelihood_bd_in_past(phylo,tot_time,time_stop, f.lamb.par,f.mu.par,desc, tot_desc,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt,cond=cond)
      
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    
    lamb.par <- temp$par[1:length(lamb_par)]
    mu.par <- temp$par[(1+length(lamb_par)):length(init)]
    f.lamb.par <- function(t){abs(f.lamb(t, lamb.par))}
    f.mu.par <- function(t){abs(f.mu(t, mu.par))}
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1) , lamb_par=lamb.par, mu_par=mu.par, f.lamb=Vectorize(f.lamb.par), f.mu=Vectorize(f.mu.par))
    
  }else{
    
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init){
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
      LH <- likelihood_bd_in_past(phylo,tot_time,time_stop, f.lamb.par,f.mu.par,desc, tot_desc,cst.lamb=cst.lamb,cst.mu=TRUE,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    f.lamb.par <- function(t){abs(f.lamb(t, lamb.par))}
    f.mu.par <- function(t){abs(f.mu(t, mu_par))}
    res <- list(model = "pure birth", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=lamb.par, f.lamb=Vectorize(f.lamb.par))
  }
  class(res) <- "fit.bd"
  return(res)
}



fit_env_in_past <- function (phylo, env_data, tot_time, time_stop, f.lamb, f.mu, lamb_par, mu_par, 
							 desc, tot_desc, df=NULL, meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE, 
							 expo.lamb=FALSE, expo.mu=FALSE, fix.mu=FALSE, dt=0, cond="crown"){
  if (tot_time > max(env_data[,1])){
    stop("The environmental data does not cover the time span of the phylogeny: either enter data that covers the full time span or run analyses on younger clades")
  }
  
  # first a spline is used to build the approximation model Env(t)
  if (is.null(df)){
    df <- smooth.spline(x=env_data[,1], env_data[,2])$df
  }
  spline_result <- sm.spline(env_data[,1],env_data[,2], df=df)
  env_func <- function(t){predict(spline_result,t)}
  # In order to perform computation, the env_func is tabulated
  # control from lower_bound -10%, upper_bound + 10%
  lower_bound_control <- 0.10
  upper_bound_control <- 0.10
  lower_bound <- min(env_data[,1])
  upper_bound <- max(env_data[,1])
  # Tabulation of the function from lower_bound -10%, upper_bound + 10%
  time_tabulated <- seq(from=lower_bound*(1.0-lower_bound_control),
                        to=upper_bound*(1.0+upper_bound_control),
                        length.out=1+1e6)
  env_tabulated <- env_func(time_tabulated)
  # Tabulated function
  env_func_tab <- function(t){
    b <- upper_bound * (1.0 + upper_bound_control)
    a <- lower_bound * (1.0 - lower_bound_control)
    # number of intervals
    n <- length(env_tabulated) - 1
    index <- 1 + as.integer( (t - a) * n / (b - a))
    return(env_tabulated[index])
  }
  f.lamb.env <- function(t,y){f.lamb(t, env_func_tab(t), y)}
  f.mu.env <- function(t,y){f.mu(t, env_func_tab(t), y)}
  res <- fit_bd_in_past(phylo, tot_time=tot_time, time_stop=time_stop, f.lamb=f.lamb.env, f.mu=f.mu.env, desc=desc, tot_desc=tot_desc, lamb_par=lamb_par, mu_par=mu_par,
                        meth=meth, cst.lamb=cst.lamb, cst.mu=cst.mu, expo.lamb=expo.lamb, expo.mu=expo.mu, fix.mu=fix.mu, dt=dt, cond=cond)
  res$model <- "environmental birth death"
  res$f.lamb <- function(t){ abs(f.lamb(t, env_func_tab(t), res$lamb_par))}
  if (fix.mu==FALSE){
    res$f.mu <- function(t){ abs(f.mu(t, env_func_tab(t), res$mu_par))}
  }
  class(res) <- "fit.env"
  return(res)
}