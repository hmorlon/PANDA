fit_bd_backbone <- function(phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = 1, 
                            backbone, spec_times, branch_times, # options for backbone analysis
                            meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, 
                            expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, dt = 1e-3, 
                            cond = "crown", model) 
{
  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"")
  
  rate <- function(par, tot_time, model){
    
    agei <- tot_time[[1]]
    
    if (grepl("BCST", model)){
      rate_data <- data.frame(Lambda = par)
      rate_spec <- rep(rate_data$Lambda[1], length(c(seq(agei,0,by=-1),0)))
      rate_ext <- rep(NA, length(c(seq(agei,0,by=-1),0)))
    }
    
    if (grepl("BCST_DCST", model)){
      rate_data <- data.frame(Lambda = par[1], Mu = par[2])
      rate_spec <- rep(rate_data$Lambda[1], length(c(seq(agei,0,by=-1),0)))
      rate_ext <- rep(rate_data$Mu[1], length(c(seq(agei,0,by=-1),0)))
    }
    
    if (grepl("BVAR", model)){
      rate_data <- data.frame(Lambda = par[1], Alpha = par[2])
      rate_spec <- rate_data$Lambda[1] *exp(rate_data$Alpha[1] * c(seq(agei,0,by=-1),0))
      rate_ext <- rep(NA, length(c(seq(agei,0,by=-1),0)))
    }
    
    if (grepl("BVAR_DCST", model)){
      rate_data <- data.frame(Lambda = par[1], Alpha = par[2], Mu = par[3])
      rate_spec <- rate_data$Lambda[1] *exp(rate_data$Alpha[1] * c(seq(agei,0,by=-1),0))
      rate_ext <- rep(rate_data$Mu[1], length(c(seq(agei,0,by=-1),0)))
    }
    
    if (grepl("BCST_DVAR", model)){
      rate_data <- data.frame(Lambda = par[1], Mu = par[2], Beta = par[3])
      rate_spec <- rep(rate_data$Lambda[1], length(c(seq(agei,0,by=-1),0)))
      rate_ext <- rate_data$Mu[1] *exp(rate_data$Beta[1] * c(seq(agei,0,by=-1),0))
    }
    
    if (grepl("BVAR_DVAR", model)){
      rate_data <- data.frame(Lambda = par[1], Alpha = par[2], Mu = par[3], Beta = par[4])
      rate_spec <- rate_data$Lambda[1] *exp(rate_data$Alpha[1] * c(seq(agei,0,by=-1),0))
      rate_ext <- rate_data$Mu[1] *exp(rate_data$Beta[1] * c(seq(agei,0,by=-1),0))
    }
    rate_df2 <- rbind(rate_spec, rate_ext)
    
    return(rate_df2)
  }
  
  if (!inherits(phylo, "phylo")) 
    stop("object \"phylo\" is not of class \"phylo\"")
  nobs <- Ntip(phylo)
  if (fix.mu == FALSE) {
    init <- c(lamb_par, mu_par)
    p <- length(init)
    optimLH <- function(init, tot_time. = tot_time, model. = model) {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1 + length(lamb_par)):length(init)]
      f.lamb.par <- function(t) {
        abs(f.lamb(t, lamb_par))
      }
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
      rate.test <- rate(init, tot_time, model) 
      
      check1 <- length(lamb_par) == 1 & length(mu_par) == 1 & lamb_par[1] <= mu_par[1]
      check2 <- min(rate.test, na.rm = T) < 0
      if(all(is.na(rate.test[2,]))){
        check3 <- F
      } else {
        check3 <- ifelse(apply(rate.test, 2, function(x) x[1] - x[2])[1] < 0, T, F)
      }
      
      if(any(c(check1, check2, check3) == T)){
        LH <- Inf
      } else {
        LH <- likelihood_bd_backbone(phylo, tot_time, f, f.lamb.par, f.mu.par,
                                   backbone, spec_times, branch_times,
                                   cst.lamb = cst.lamb, cst.mu = cst.mu, 
                                   expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt, 
                                   cond = cond)
      }
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    mu.par <- temp$par[(1 + length(lamb_par)):length(init)]
    f.lamb.par <- function(t) {
      abs(f.lamb(t, lamb.par))
    }
    f.mu.par <- function(t) {
      abs(f.mu(t, mu.par))
    }
    res <- list(model = "birth death", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1))/(nobs - 
                                                                     p - 1), lamb_par = lamb.par, mu_par = mu.par, 
                f.lamb = Vectorize(f.lamb.par), f.mu = Vectorize(f.mu.par))
  }
  else {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init, tot_time. = tot_time, model. = model) {
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t) {
        abs(f.lamb(t, lamb_par))
      }
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
      
      rate.test <- rate(init, tot_time, model) 
      check2 <- min(rate.test, na.rm = T) < 0
      
      if(check2){
        LH <- Inf
      } else {
        
      LH <- likelihood_bd_backbone(phylo, tot_time, f, f.lamb.par, f.mu.par, 
                                   backbone, spec_times, branch_times,
                                   cst.lamb = cst.lamb, cst.mu = TRUE, 
                                   expo.lamb = expo.lamb, dt = dt, cond = cond)
      }
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    f.lamb.par <- function(t) {
      abs(f.lamb(t, lamb.par))
    }
    f.mu.par <- function(t) {
      abs(f.mu(t, mu_par))
    }
    res <- list(model = "birth.death", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1))/(nobs - 
                                                                     p - 1), lamb_par = lamb.par, f.lamb = Vectorize(f.lamb.par))
  }
  class(res) <- "fit.bd"
  return(res)
}
