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
nameEx("InfTemp")
### * InfTemp

flush(stderr()); flush(stdout())

### Name: InfTemp
### Title: Inferred Temperature
### Aliases: InfTemp
### Keywords: datasets

### ** Examples

data(InfTemp)
plot(InfTemp)



cleanEx()
nameEx("Phocoenidae")
### * Phocoenidae

flush(stderr()); flush(stdout())

### Name: Phocoenidae
### Title: Phocoenidae phylogeny
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
# Fit of a pure birth model (no extinction) with an
# exponential variation of speciation rate with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
result_exp <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89,
                     expo.lamb = TRUE, fix.mu=TRUE)
result_exp$model <- "birth death - exponential variation rate"

# Fit of a pure birth model (no extinction) with a
# linear variation of speciation rate with time
f.lamb <-function(t,y){y[1] + y[2] * t}
f.mu<-function(t,y){0}
lamb_par<-c(0.09, 0.001)
mu_par<-c()
result_lin <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, fix.mu=TRUE)
result_lin$model <- "birth death - linear variation rate"

# Fit of a pure birth model (no extinction) with a constant speciation rate
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.09, 0.001)
mu_par<-c()
result_cst <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89,
                     cst.lamb = TRUE, fix.mu=TRUE)
result_cst$model <- "birth death - cst rate"

# Find the best model
index <- which.min(c(result_exp$aicc, result_lin$aicc, result_cst$aicc))
rbind(result_exp, result_lin, result_cst)[index,]



cleanEx()
nameEx("fit_coal_cst")
### * fit_coal_cst

flush(stderr()); flush(stdout())

### Name: fit_coal_cst
### Title: Maximum likelihood fit of the equilibrium model
### Aliases: fit_coal_cst

### ** Examples

data(Cetacea)
result <- fit_coal_cst(Cetacea, tau0=1.e-3, gamma=-1, cst.rate=FALSE, N0=89)
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
result <- fit_coal_var(Cetacea, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0, N0=89)
print(result)



cleanEx()
nameEx("fit_env_bd")
### * fit_env_bd

flush(stderr()); flush(stdout())

### Name: fit_env_bd
### Title: Maximum likelihood fit of the environmental birth-death model
### Aliases: fit_env_bd

### ** Examples

data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
data(InfTemp)
df<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

# fits a model with lambda varying as an exponential function of temperature and mu fixed to 0 (no extinction)

f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
result_exp <- fit_env_bd(Cetacea,InfTemp,tot_time,f.lamb,f.mu,
      lamb_par,mu_par,f=87/89, fix.mu=TRUE, df=df)



cleanEx()
nameEx("likelihood_bd")
### * likelihood_bd

flush(stderr()); flush(stdout())

### Name: likelihood_bd
### Title: Likelihood of a phylogeny under the general birth-death model
### Aliases: likelihood_bd

### ** Examples

data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
# Compute the likelihood for a pure birth model (no extinction) with
# an exponential variation of speciation rate with time
lamb_par <- c(0.1, 0.01)
f.lamb <- function(t){lamb_par[1] * exp(lamb_par[2] * t)}
f.mu <- function(t){0}
f <- 87/89
likelihood_bd(Cetacea,tot_time,f.lamb,f.mu,f,cst.mu=TRUE,expo.lamb=TRUE)



cleanEx()
nameEx("likelihood_coal_cst")
### * likelihood_coal_cst

flush(stderr()); flush(stdout())

### Name: likelihood_coal_cst
### Title: Likelihood of a phylogeny under the equilibrium diversity model
### Aliases: likelihood_coal_cst

### ** Examples

data(Cetacea)
Vtimes <- sort(branching.times(Cetacea))
tau0 <- 0.1
gamma <- 0.001
ntips <- Ntip(Cetacea)
N0 <- 89
likelihood <- likelihood_coal_cst(Vtimes,ntips,tau0, gamma,N0)



cleanEx()
nameEx("likelihood_coal_var")
### * likelihood_coal_var

flush(stderr()); flush(stdout())

### Name: likelihood_coal_var
### Title: Likelihood of a birth-death model using a coalescent approch
### Aliases: likelihood_coal_var

### ** Examples

data(Cetacea)
Vtimes <- sort(branching.times(Cetacea))
lamb0 <- 0.1
alpha <- 0.001
mu0<-0
beta<-0
ntips <- Ntip(Cetacea)
N0 <- 89
likelihood <- likelihood_coal_var(Vtimes, ntips, lamb0, alpha, mu0, beta, N0)



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
