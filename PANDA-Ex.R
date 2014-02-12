pkgname <- "PANDA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PANDA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Cetacea")
### * Cetacea

flush(stderr()); flush(stdout())

### Name: Cetacea
### Title: Cetacean phylogeny
### Aliases: Cetacea
### Keywords: datasets

### ** Examples

data(Cetacea)
print(Cetacea)
plot(Cetacea)



cleanEx()
nameEx("Phocoenidae")
### * Phocoenidae

flush(stderr()); flush(stdout())

### Name: Phocoenidae
### Title: Phocoenidae
### Aliases: Phocoenidae
### Keywords: datasets

### ** Examples

data(Phocoenidae)
print(Phocoenidae)
plot(Phocoenidae)



cleanEx()
nameEx("fit_bd")
### * fit_bd

flush(stderr()); flush(stdout())

### Name: fit_bd
### Title: Maximum likelihood fit of the general birth-death model
### Aliases: fit_bd

### ** Examples

data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
#fit of a pure birth model (no extinction) with an exponential variation of speciation rate with time
f.lamb <-function(x,y){y[1] * exp(y[2] * x)}
f.mu<-function(x,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, expo.lamb = TRUE, fix.mu=TRUE)



cleanEx()
nameEx("fit_coal_cst_mod")
### * fit_coal_cst_mod

flush(stderr()); flush(stdout())

### Name: fit_coal_cst_mod
### Title: Fit diversity model with time-constant rate
### Aliases: fit_coal_cst_mod

### ** Examples

data(Cetacea)
result <- fit_coal_cst_mod(Cetacea, tau0=1.5, N0=1e6)
print(result)



cleanEx()
nameEx("fit_coal_exp_mod")
### * fit_coal_exp_mod

flush(stderr()); flush(stdout())

### Name: fit_coal_exp_mod
### Title: Fit diversity model with time-varying rate
### Aliases: fit_coal_exp_mod

### ** Examples

data(Cetacea)
result <- fit_coal_exp_mod(Cetacea, tau0=1.e-3, gamma=-1)
print(result)



cleanEx()
nameEx("fit_coal_var")
### * fit_coal_var

flush(stderr()); flush(stdout())

### Name: fit_coal_var
### Title: Fit birth-death model using a coalescent approch
### Aliases: fit_coal_var

### ** Examples

data(Cetacea)
result <- fit_coal_var(Cetacea, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0, N0=1e6)
print(result)



cleanEx()
nameEx("likelihood_bd")
### * likelihood_bd

flush(stderr()); flush(stdout())

### Name: likelihood_bd
### Title: Compute the likelihood of a phylogeny under the general
###   birth-death model
### Aliases: likelihood_bd

### ** Examples

data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
#compute the likelihood for a pure birth model (no extinction) with an exponential variation of speciation rate with time
lamb_par <- c(0.1, 0.01)
f.lamb <- function(x){lamb_par[1] * exp(lamb_par[2] * x)}
f.mu <- function(x){0}
f <- 87/89
likelihood_bd(Cetacea,tot_time,f.lamb,f.mu,f,cst.mu=TRUE,expo.lamb=TRUE)



cleanEx()
nameEx("likelihood_coal_cst_mod")
### * likelihood_coal_cst_mod

flush(stderr()); flush(stdout())

### Name: likelihood_coal_cst_mod
### Title: Likelihood of a constant diversity model with time-constant rate
### Aliases: likelihood_coal_cst_mod

### ** Examples

data(bird.families)
Vtimes <- sort(branching.times(bird.families))
tau0 <- 0.1
ntips <- Ntip(bird.families)
N0 <- 1e6
likelihood <- likelihood_coal_cst_mod(Vtimes,ntips,tau0,N0)



cleanEx()
nameEx("likelihood_coal_exp_mod")
### * likelihood_coal_exp_mod

flush(stderr()); flush(stdout())

### Name: likelihood_coal_exp_mod
### Title: Likelihood of a constant diversity model with time-varying rate
### Aliases: likelihood_coal_exp_mod

### ** Examples

data(bird.families)
Vtimes <- sort(branching.times(bird.families))
tau0 <- 0.1
gamma <- -0.001
ntips <- Ntip(bird.families)
N0 <- 1e6
likelihood <- likelihood_coal_exp_mod(Vtimes,ntips,tau0, gamma,N0)



cleanEx()
nameEx("likelihood_coal_var")
### * likelihood_coal_var

flush(stderr()); flush(stdout())

### Name: likelihood_coal_var
### Title: Likelihood of a birth-death model using a coalescent approch
### Aliases: likelihood_coal_var

### ** Examples

data(bird.families)
Vtimes <- sort(branching.times(bird.families))
# likelihood considering model 3
lamb0 <- 0.1
alpha <- -0.001
ntips <- Ntip(bird.families)
likelihood <- likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0=0, beta=0, N0=1e6)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
