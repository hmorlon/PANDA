## ------------
## Proposal kernels
## ------------

proposal_kernel <- function(value, tuning, method){
	
	n <- length(value)
	hastings=0
	# 3.4641016 is sqrt(12) see Yang & Rodriguez 2013
	switch(method,
	"bactrian"={
		newval = value + rbactrian(n)*tuning
	},
    "bactrianTriangle"={
        newval = value + rbactrianTriangle(n)*tuning
    },
    "bactrianLaplace"={
        newval = value + rbactrianLaplace(n)*tuning
    },
	"uniform"={
		newval = value + (runif(n,0,1)-0.5)*3.4641016*tuning;
	},
    "uniform-sliding"={
        newval = value + (runif(n,0,1)-0.5)*3.4641016*tuning;
        newval = reflect(newval,0,100)
    },
    "sliding"={
        newval = sliding_rate(value, tuning=tuning, min=0, max=Inf)
    },
	"Muniform"={
	    hastings =  (runif(n,0,1)-0.5)*3.4641016*tuning
        newval = value * exp(hastings);
        # then sum the hastings ratio because we use n multiplier on the log-scale
        hastings = sum(hastings)
	},
    "MuniformLog"={
        hastings =  log(value)+(runif(n,0,1)-0.5)*3.4641016*tuning
        hastings = reflect(hastings, -99, 99)
        newval = exp(hastings);
        # then sum the hastings ratio because we use n multiplier on the log-scale
        hastings = sum(hastings - log(value))
    },
	"Mbactrian"={
		hastings =  rbactrian(n)*tuning
        newval = value * exp(hastings);
        # hastings ratio
        hastings = sum(hastings)
	},
    "MbactrianTriangle"={
        hastings =  rbactrianTriangle(n)*tuning
        newval = value * exp(hastings);
        # hastings ratio
        hastings = sum(hastings)
    },
    "MbactrianLaplace"={
        hastings =  rbactrianLaplace(n)*tuning
        newval = value * exp(hastings);
        # hastings ratio
        hastings = sum(hastings)
    },
    "MbactrianLog"={
        hastings =  log(value) + rbactrian(n)*tuning
        hastings = reflect(hastings, -99, 99)
        newval = exp(hastings);
        # must add bounds reflexions here (e.g. between -99 and 99)
        # hastings ratio
        hastings = sum(hastings - log(value))
    },
    "multiplier"={
        prop = mult_rate(value, tuning)
        newval = prop[1]
        hastings = prop[2]
    },

	stop("Proposal kernel not yet implemented!")
	)
	
	results <- list(moves=newval, hastings=hastings)
	return(results)
}


## ------- Multiplier proposal of Lakner et al. 2008 - Systematic Biology
mult_rate <-
function( rate, tuning=1.5){
  tuning = 2*log(tuning) # m will be in the interval [1/tu.; tu.]
  u = runif(1,0,1)
  h = tuning*(u-0.5) # the Hastings ratio
  m = exp(h)
  rate = rate * m
  
  return(c(rate, h))
}


## -------- Sliding windows proposal with reflexion bounds
sliding_rate <-
function( rate, tuning=1.5, min=0, max=Inf){
    n <- length(rate)
  repeat{
    u=runif(n,0,1)
    v=rate+(u-0.5)*3.4641016*tuning # See Yang & Rodriguez 2013
    
    if(any(v>max)) {
      v[v>max]=max-(v[v>max]-max)
    }
    if(any(v<min)){
      v[v<min]=min-(v[v<min]-min)
    }
    # break the loop
    if(any(v<min)==FALSE & any(v>max)==FALSE) break;
  }
  # Hastings ratio is 0	?;
  return(v)
}


## --------- Reflection for bounded parameters from Yang & Rodriguez 2013
reflect <- function(param, min_val, max_val){
    e=0
    if(param < min_val) {e = min_val - param; side = 0}
    if(param > max_val) {e = param - max_val; side = 1}
    if(e){
        n = as.integer(e/(max_val - min_val))
        if(n%%2 == 1) side = 1-side
        e = e - n*(max_val - min_val)
        param = ifelse(side, max_val-e, min_val+e)
    }
    return(param)
}



# Geometric Brownian motion expectation - eq. 10 in Guindon 2013 - Systematic Biology
geometricBM <- function(par, parents, times, sigma, len){
  results <- .Call("geometricExpectation", lambda=par, parents=as.integer(parents), length=as.integer(len), sigma=sigma, brlength=times)
  return(results)
}

# Geometric Brownian motion from a Gamma distribution parameterized by the moments - eq. 10, 24 in Guindon 2013 - Systematic Biology
geometricBMGibbs <- function(par, parents, times, sigma, len){
  results <- .Call("geometricExpectationGibbs", lambda=par, parents=as.integer(parents), length=as.integer(len), sigma=sigma, brlength=times)
  return(results)
}

# Arithmetic average of nodes-rates; e.g. Kishino et al. 2001
arithmeticBM <- function(par, parents, len){
  results <- .Call("arithmetic", lambda=par, parents=as.integer(parents), length=as.integer(len))
  return(results)
}

# Simulate a normal bactrian variate (Yang & Rodriguez 2013 - PNAS)

rbactrian <- function(n, m=0.95){
  mBactrian = m
  sBactrian = sqrt(1-m^2)
  z = mBactrian + rnorm(n,0,1)*sBactrian
  rdunif <- runif(n,0,1)<0.5
  sign <- ifelse(rdunif,-1,1)
  z=z*sign
  return(z)
}

# Simulate a triangular bactrian variate (Yang & Rodriguez 2013 - PNAS)

rbactrianTriangle <- function(n, m=0.95){
  mBactrian = m
  sBactrian = sqrt(1-m^2)
  # triangle variate
  u <- runif(n,0,1)
  usign <- u<0.5
  variate = ifelse(usign,sqrt(6)-2*sqrt(3*(1-u)),-sqrt(6)+2*sqrt(3*u))
  
  z = mBactrian + variate*sBactrian
  rdunif <- runif(n,0,1)<0.5
  sign <- ifelse(rdunif,-1,1)
  z=z*sign
  return(z)
}

# Simulate a laplace bactrian variate (Yang & Rodriguez 2013 - PNAS)

rbactrianLaplace <- function(n, m=0.95){
  mBactrian = m
  sBactrian = sqrt(1-m^2)
  u <- runif(n,0,1) - 0.5
  variate = log(1-2*abs(u)) * 0.70710678118654752440
  variate = ifelse(u>=0, -variate, variate)
  
  z = mBactrian + variate*sBactrian
  rdunif <- runif(n,0,1)<0.5
  sign <- ifelse(rdunif,-1,1)
  z=z*sign
  return(z)
}
