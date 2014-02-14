.Phi <- function(t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE)
{

  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- 1-exp(r*t)/(1/f+lamb/r*(exp(r*t)-1))
    return(res)
  }

  if ((cst.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }

  if ((expo.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }

  if ((expo.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0/beta*(exp(beta*y)-exp(beta*x))}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }

  else
  {
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int <- function(x,y){.Integrate(Vectorize(r),x,y,stop.on.error=FALSE)}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    res
  }
}
