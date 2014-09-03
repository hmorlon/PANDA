plot.fit.bd <- function(fit.bd,tot_time)
{
  t <- seq(0,tot_time, length.out=100)
  X11()
  plot(t, fit.bd$f.lamb(t), type='l', xlab="Time", ylab="Speciation rate", xlim=c(tot_time,0),main="Fitted speciation rate")
  X11()
  plot(t, fit.bd$f.mu(t), type='l', xlab="Time", ylab="Extinction rate", xlim=c(tot_time,0),main="Fitted extinction rate")
  r <- function(t) {fit.bd$f.lamb(t) - fit.bd$f.mu(t)}
  X11()
  plot(t, r(t), type='l', xlab="Time", ylab="Diversification rate", xlim=c(tot_time,0),main="Fitted diversification rate")
}
