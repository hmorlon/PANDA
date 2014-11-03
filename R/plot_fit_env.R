plot_fit <- function(x, env_data, tot_time)
{
  if (!inherits(x, "fit.env.bd"))
      stop("object is not of class \"fit.env.bd\"")
  t <- seq(0,tot_time, length.out=100)
  dev.new()
  plot(-t, x$f.lamb(t), type='l', xlab="time", ylab="speciation rate", main="Fitted speciation rate")
  # Plot f.lamb(env_data)
  df <- smooth.spline(x=env_data[,1], env_data[,2])$df
  spline_result <- sm.spline(env_data[,1],env_data[,2], df=df)
  env_func <- function(t){predict(spline_result,t)}
  dev.new()
  plot(env_func(t), x$f.lamb(t), type='l', xlab="Environmental data", ylab="speciation rate", main="Fitted speciation rate")


  if ("f.mu" %in% attributes(x))
  {
    # Attribute f.mu ==> not fixed extinction
    dev.new()
    plot(-t, x$f.mu(t), type='l', xlab="time", ylab="extinction rate", main="Fitted extinction rate")
    plot(env_func(t), x$f.mu(t), type='l', xlab="Environmental data", ylab="extinction rate", main="Fitted extinction rate")
    r <- function(t) {x$f.lamb(t) - x$f.mu(t)}
    dev.new()
    plot(-t, r(t), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
    plot(env_func(t), r(t), type='l', xlab="Environmental data", ylab="net diversification rate", main="Fitted net diversification rate")
  }
  else
  {
    dev.new()
    plot(-t, x$f.lamb(t), type='l', xlab="time", ylab="net diversification rate", main="Fitted net diversification rate")
    plot(env_func(t), x$f.lamb(t), type='l', xlab="Environmental data", ylab="net diversification rate", main="Fitted net diversification rate")
  }
}
