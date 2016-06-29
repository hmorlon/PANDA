fit_bd <-
 function (phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1,
           meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE,
           expo.lamb=FALSE, expo.mu=FALSE, fix.mu=FALSE,
           dt=0, cond="crown", mcmc = F, mcmcSettings = NULL, prior = NULL)
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  nobs <- Ntip(phylo)

  if (fix.mu==FALSE)
  {
    init <- c(lamb_par,mu_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1+length(lamb_par)):length(init)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt=dt,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    mu.par <- temp$par[(1+length(lamb_par)):length(init)]
    f.lamb.par <- function(t){f.lamb(t, lamb.par)}
    f.mu.par <- function(t){f.mu(t, mu.par)}
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1) , lamb_par=lamb.par, mu_par=mu.par, f.lamb=Vectorize(f.lamb.par), f.mu=Vectorize(f.mu.par))
  }

  else
  {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=TRUE,expo.lamb=expo.lamb,dt=dt,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    f.lamb.par <- function(t){f.lamb(t, lamb.par)}
    f.mu.par <- function(t){f.mu(t, mu_par)}
    res <- list(model = "birth.death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=lamb.par, f.lamb=Vectorize(f.lamb.par))
  }
  
  if (mcmc==TRUE){
    
    require(BayesianTools) # Note Florian: I attach the BT package at this point to make sure that the BT analysis functions are available to the user after the MCMCs have been run. We could get rid of this of course, but I think it's convenient has has minimal risk, lower at any rate than attaching BT with the package. If you would remove this, users would have to load library BT to analyze the MCMC output. 
    
    if(fix.mu==FALSE) mcmcSettings$startValue = c(lamb.par, mu.par)
    LL <- function(par) -optimLH(par)
    names = paste("lambda", 1:length(lamb.par))
    if (fix.mu==FALSE) names = c(names, paste("mu", 1:length(mu.par)))
    if(is.null(prior)) prior = BayesianTools::createUniformPrior(lower = rep(0,p), upper = rep(20,p))
    BayesianSetup = BayesianTools::createBayesianSetup(likelihood = LL, prior = prior, names = names)
    MCMCres<-BayesianTools::runMCMC(bayesianSetup = BayesianSetup, sampler = "DEzs", settings = mcmcSettings)
    res$mcmc = MCMCres
  }
  
  class(res) <- "fit.bd"
  return(res)
}
