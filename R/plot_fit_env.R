plot_fit_env <- function(fit.env, env_data, tot_time)
{
  if (!inherits(fit.env, "fit.env"))
      stop("object is not of class \"fit.env\"")
  t <- seq(0,tot_time, length.out=100)
  dev.new()
  plot(-t, sapply(t,fit.env$f.lamb), type='l', xlab="time", ylab="speciation rate", main="Fitted speciation rate")
  # Plot f.lamb(env_data)
  df <- smooth.spline(env_data[,1], env_data[,2])$df
  spline_result <- sm.spline(env_data[,1],env_data[,2], df=df)
  env_func <- function(t){predict(spline_result,t)}
  dev.new()
  plot(env_func(t), sapply(t,fit.env$f.lamb), xlab="Environmental data", ylab="speciation rate", main="Fitted speciation rate")


  if ("f.mu" %in% attributes(fit.env)$names)
  {
    # Attribute f.mu ==> not fixed extinction
    dev.new()
    plot(-t, sapply(t,fit.env$f.mu), type='l', xlab="time", ylab="extinction rate", main="Fitted extinction rate")
    dev.new()
    plot(env_func(t), sapply(t,fit.env$f.mu), xlab="Environmental data", ylab="extinction rate", main="Fitted extinction rate")
    r <- function(t) {fit.env$f.lamb(t) - fit.env$f.mu(t)}
    dev.new()
    plot(-t, sapply(t,r), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
  dev.new()
    plot(env_func(t), sapply(t,r), xlab="Environmental data", ylab="net diversification rate", main="Fitted net diversification rate")
  }
  else
  {
    dev.new()
    plot(-t, sapply(t,fit.env$f.lamb), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
    dev.new()
    plot(env_func(t), sapply(t,fit.env$f.lamb), xlab="Environmental data", ylab="net diversification rate", main="Fitted net diversification rate")
  }
}
