# ---- Function to simulate random covariance matrices with a specified eigenstructure
# From Uyeda et al. 2015 - Systematic Biology 64(4):677-689.
Posdef <- function (p, ev = rexp(p, 1/100)) {
  Z <- matrix(ncol=p, rnorm(p^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}
