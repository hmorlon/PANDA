.defineRintRintint <- function(f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=1e-3,tot_time)
{
  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    r.int <- function(x,y){r*(y-x)}
    r.int.int <- function(x,y){lamb*(exp(r*y)-exp(r*x))/r}
    return(list(r.int=r.int, r.int.int=r.int.int))
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
    return(list(r.int=r.int, r.int.int=r.int.int))
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
    return(list(r.int=r.int, r.int.int=r.int.int))
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
    return(list(r.int=r.int, r.int.int=r.int.int))
  }

  else
  {
    ageMin <- 0
    ageMax <- tot_time
    # Constant approximation on steps of size dt
    Nintervals <- 1 + as.integer((ageMax-ageMin)/dt)
    X <- seq(ageMin, ageMax, length.out = Nintervals + 1)
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int.tab <- cumsum(r(X)) * (ageMax - ageMin) / Nintervals
    r.int <- function(x,y)
    {
      indy <- 1 + as.integer( (y - ageMin) * Nintervals / (ageMax - ageMin))
      indx <- 1 + as.integer( (x - ageMin) * Nintervals / (ageMax - ageMin))
      value <- r.int.tab[indy] - r.int.tab[indx]
      return(value)
    }
    r.int.0 <- function(y){exp(r.int.tab[1 + as.integer( (y - ageMin) * Nintervals / (ageMax - ageMin))]) * f.lamb(y)}
    r.int.int.tab <- cumsum(r.int.0(X)) * (ageMax - ageMin) / Nintervals
    r.int.int <- function(x,y)
    {
      indy <- 1 + as.integer( (y - ageMin) * Nintervals / (ageMax - ageMin))
      indx <- 1 + as.integer( (x - ageMin) * Nintervals / (ageMax - ageMin))
      value <- r.int.int.tab[indy] - r.int.int.tab[indx]
      return(value)
    }
    return(list(r.int=r.int, r.int.int=r.int.int))
  }
}
