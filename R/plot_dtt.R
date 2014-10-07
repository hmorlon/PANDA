plot_dtt <- function(fit.bd,tot_time,N0)
{
  t <- seq(0,tot_time, length.out=100)
  r <- function(t) {-fit.bd$f.lamb(t) + fit.bd$f.mu(t)}
  R <- function(s){.Integrate(r,0,s)}
  N <- N0 * exp(Vectorize(R)(t))
  dev.new()
  plot(t, N, type='l', xlab="time", ylab="number of species", xlim=c(tot_time,0),main="Diversity Through Time")
}
