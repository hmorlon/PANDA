fit_env_bd <-
 function (phylo, env_data, tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1,
           meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE,
           expo.lamb=FALSE, expo.mu=FALSE, fix.mu=FALSE,
           cond="crown")
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  nobs <- Ntip(phylo)

  # first a spline is used to build the approximation model Env(t)
  spline_result <- sm.spline(env_data[,1],env_data[,2])
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
                        length.out=floor(1.e6))
  env_tabulated <- env_func(time_tabulated)
  # Tabulated function
  env_func_tab <- function(t)
  {
    b <- upper_bound * (1.0 + upper_bound_control)
    a <- lower_bound * (1.0 - lower_bound_control)
    n <- length(env_tabulated)
    index <- 1 + as.integer( (t - a) * n / (b - a))
    return(env_tabulated[index])
  }
  if (fix.mu==FALSE)
  {
    init <- c(lamb_par,mu_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1+length(lamb_par)):length(init)]
      # function lambda & mu must take into account env_data
      # We redefine the function f.lamb & f.mu
      f.lamb.par <- function(t){abs(f.lamb(env_func_tab(t),lamb_par))}
      f.mu.par <- function(t){abs(f.mu(env_func_tab(t),mu_par))}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,cond=cond)
      print(c("LH",LH))
      return(-LH)
    }
    temp <- optim(init, optimLH, method = meth)
    res <- list(model = "environmental birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1) , lamb_par=temp$par[1:length(lamb_par)], mu_par=temp$par[(1+length(lamb_par)):length(init)])
  }

  else
  {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      # function lambda & mu must take into account env_data
      # We redefine the function f.lamb & f.mu
      f.lamb.par <- function(t){abs(f.lamb(env_func_tab(t),lamb_par))}
      f.mu.par <- function(t){abs(f.mu(env_func_tab(t),mu_par))}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=TRUE,expo.lamb=expo.lamb,expo.mu=FALSE,cond=cond)
      return(-LH)
    }
    temp <- optim(init, optimLH, method = meth)
    res <- list(model = "environmental birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$par[1:length(lamb_par)])
  }
  return(res)
}
