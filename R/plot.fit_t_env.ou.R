################################################################################
##                                                                            ##
##                         RPANDA : Climatic-model                            ##
##                                                                            ##
##   Julien Clavel - 25-01-2024                                               ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################

plot.fit_t.env.ou<-function(x, steps=100, ...){
  
  # Rates through time function
  alpha = x$param$alpha
  param = x$param$par
  theta = x$param$theta_0
  fun_temp = x$env_func
  model = x$model
  root_age = x$tot_time
  
  # Times steps
  time <- seq(0, root_age, length.out=steps)
  
  # Rates through time
  fun_ou_expectation <- function(x, stop_time, par, theta_0, fun) { alpha*exp(alpha*(x-stop_time))*model(root_age - x, fun, par, theta_0) } 
  # fun_temp is the function used for model fit
  
  ou_clim_expectation = sapply(time, function(ti) exp(-alpha*ti)*model(root_age, env=fun_temp, param=param, theta0=theta) + integrate(fun_ou_expectation, 0, ti, stop_time=ti, par=param, theta_0=theta, fun=fun_temp, subdivisions=200L, rel.tol = .Machine$double.eps^0.05)$value) 
  # alpha; beta; theta are the parameters estimated by the model
  
  # plot the curves
  plot(-time, rev(ou_clim_expectation), type="l", las=1, xlab="Times", ylab=bquote(paste("Evolutionary trajectory ", theta)), ...)
  
  results<-list(time_steps=time, values=ou_clim_expectation)
  invisible(results)
}

# Allows drawing lines and superposing various results

lines.fit_t.env.ou<-function(x,steps=100,...){
  
  # Rates through time function
  alpha = x$param$alpha
  param = x$param$par
  theta = x$param$theta_0
  fun_temp = x$env_func
  model = x$model
  root_age = x$tot_time
  
  # Times steps
  time <- seq(0, root_age, length.out=steps)
  
  # Rates through time
  fun_ou_expectation <- function(x, stop_time, par, theta_0, fun) { alpha*exp(alpha*(x-stop_time))*model(root_age - x, fun, par, theta_0) } 
  # fun_temp is the function used for model fit
  
  ou_clim_expectation = sapply(time, function(ti) exp(-alpha*ti)*model(root_age, env=fun_temp, param=param, theta0=theta) + integrate(fun_ou_expectation, 0, ti, stop_time=ti, par=param, theta_0=theta, fun=fun_temp, subdivisions=200L, rel.tol = .Machine$double.eps^0.05)$value) 
  # alpha; beta; theta are the parameters estimated by the model
  
  # plot the curves

  lines(-time, rev(ou_clim_expectation), type='l', ...)
  results<-list(time_steps=time, values=ou_clim_expectation)
  invisible(results)
}
