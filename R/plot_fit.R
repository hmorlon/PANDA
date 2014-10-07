plot_fit <- function(fit.bd,tot_time)
{
  if (!inherits(fit.bd, "fit.bd"))
      stop("object \"fit.bd\" is not of class \"fit.bd\"")
  t <- seq(0,tot_time, length.out=100)
  dev.new()
  plot(t, fit.bd$f.lamb(t), type='l', xlab="time", ylab="speciation rate", xlim=c(tot_time,0),main="Fitted speciation rate")
  dev.new()
  plot(t, fit.bd$f.mu(t), type='l', xlab="time", ylab="extinction rate", xlim=c(tot_time,0),main="Fitted extinction rate")
  r <- function(t) {fit.bd$f.lamb(t) - fit.bd$f.mu(t)}
  dev.new()
  plot(t, r(t), type='l', xlab="time", ylab="net diversification rate", xlim=c(tot_time,0),main="Fitted net diversification rate")
}
