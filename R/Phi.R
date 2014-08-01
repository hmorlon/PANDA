.Phi <- function(t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE)
{

  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- 1-r*exp(r*t)/(r/f+lamb*(exp(r*t)-1))
    return(res)
  }

  if ((cst.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }

  if ((expo.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
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
    g <- function(y){r.int(0,y)}
    gvect <- function(y){mapply(g,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- 1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
    return(res)
  }

  else
  {
    ageMin <- 0
    ageMax <- t
    dt <- 0.005
    Nintervals <- 1 + as.integer((ageMax-ageMin)/dt)
    X <- seq(ageMin, ageMax, length.out = Nintervals + 1)
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int <- cumsum(r(X)) * (ageMax - ageMin) / Nintervals
    r.int.0 <- function(y){exp(r.int[1 + as.integer( (y - ageMin) * Nintervals / (ageMax - ageMin))]) * f.lamb(y)}
    r.int.int.tab <- cumsum(r.int.0(X)) * (ageMax - ageMin) / Nintervals
    r.int.int <- function(x,y)
    {
      indy <- 1 + as.integer( (y - ageMin) * Nintervals / (ageMax - ageMin))
      indx <- 1 + as.integer( (x - ageMin) * Nintervals / (ageMax - ageMin))
      value <- r.int.int.tab[indy] - r.int.int.tab[indx]
      return(value)
    }
    rit <- r.int[1 + Nintervals]
    ri0t <- r.int.int(0,t)
    res <- 1.0 - exp(rit)/(1/f+ri0t)
    return(res)
  }
}
