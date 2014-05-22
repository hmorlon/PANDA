.Psi <- function(s,t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE)
{
  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- exp(r*(t-s))*(abs(1+(lamb*(exp(r*t)-exp(r*s)))/(r/f+lamb*(exp(r*s)-1))))^(-2)
    return(res)
  }

  ####### exponential dependencies ########

  if ((cst.lamb==TRUE) & (expo.mu==TRUE))

  {
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
    return(res)
  }

  if ((expo.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    r.int <- function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
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
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
    return(res)
  }

####### other dependencies ########

  else
  {
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int <- function(x,y){.Integrate(r,x,y,stop.on.error=FALSE)}
    r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
    r.int.int <- function(x,y){.Integrate(r.int.0,x,y,stop.on.error=FALSE)}
    rst <- r.int(s,t)
    rist <- r.int.int(s,t)
    ri0s <- r.int.int(0,s)
    if (is.infinite(rst) || is.infinite(rist) || is.infinite(ri0s))
    {
      rvect <- function(t){mapply(r,t)}
      r.int <- function(x,y){.Integrate(rvect,x,y,stop.on.error=FALSE)}
      r.int.0 <- function(y){exp(r.int(0,y))*f.lamb(y)}
      rvect0 <- function(y){mapply(r.int.0,y)}
      r.int.int <- function(x,y){.Integrate(rvect0,x,y,stop.on.error=FALSE)}
      res <- exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
      return(res)
    }
    else
    {
      res <- exp(rst)*(abs(1+rist/(1/f+ri0s)))^(-2)
      return(res)
    }
  }
}

