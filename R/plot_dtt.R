plot_dtt <- function(fit.bd,tot_time, N0=1)
{
  t <- seq(0,tot_time, length.out=100)
  r <- function(t) {-fit.bd$f.lamb(t) + fit.bd$f.mu(t)}
  R <- function(s){.Integrate(r,0,s)}
  N <- N0 * exp(Vectorize(R)(t))
  X11()
  plot(t, N, type='l', xlab="Time", ylab="dtt", xlim=c(tot_time,0),main="Diversification Through Time")
  grid()
}
