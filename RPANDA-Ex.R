pkgname <- "RPANDA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('RPANDA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Anolis.data")
### * Anolis.data

flush(stderr()); flush(stdout())

### Name: Anolis.data
### Title: Anolis dataset
### Aliases: Anolis.data

### ** Examples

data(Anolis.data)
plot(Anolis.data$phylo)
print(Anolis.data$data)
print(Anolis.data$geography.object[[16]])

# Compute the likelihood that the S value is twice the ML estimate
par <- c(0.0003139751, (2*-0.06387258))
lh <- -likelihood_MC(phylo,pPC1,par)



cleanEx()
nameEx("BICompare")
### * BICompare

flush(stderr()); flush(stdout())

### Name: BICompare
### Title: Identify modalities in a phylogeny
### Aliases: BICompare

### ** Examples

data(Cetacea)
BICompare(Cetacea,5)



cleanEx()
nameEx("Balaenopteridae")
### * Balaenopteridae

flush(stderr()); flush(stdout())

### Name: Balaenopteridae
### Title: Balaenopteridae phylogeny
### Aliases: Balaenopteridae
### Keywords: datasets

### ** Examples

data(Balaenopteridae)
print(Balaenopteridae)
plot(Balaenopteridae)



cleanEx()
nameEx("Calomys")
### * Calomys

flush(stderr()); flush(stdout())

### Name: Calomys
### Title: Calomys phylogeny
### Aliases: Calomys
### Keywords: datasets

### ** Examples

data(Calomys)
print(Calomys)
plot(Calomys)



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
nameEx("CreateGeoObject")
### * CreateGeoObject

flush(stderr()); flush(stdout())

### Name: CreateGeoObject
### Title: Create blank geography object
### Aliases: CreateGeoObject

### ** Examples


#To create a geography.object
tree<-rcoal(5)
tree_internal.nodes<-CreateGeoObject(tree)
geography.matrix<-tree_internal.nodes$geography.object
geography.matrix[[1]][,1]<-c(1,0)
geography.matrix[[1]][,2]<-c(0,1)
geography.matrix[[2]][,1]<-c(1,1,0)
geography.matrix[[2]][,2]<-c(1,1,0)
geography.matrix[[2]][,3]<-c(0,0,1)
geography.matrix[[3]][,1]<-c(1,1,0,0)
geography.matrix[[3]][,2]<-c(1,1,0,0)
geography.matrix[[3]][,3]<-c(0,0,1,1)
geography.matrix[[3]][,4]<-c(0,0,1,1)
geography.matrix[[4]][,1]<-c(1,1,1,0,0)
geography.matrix[[4]][,2]<-c(1,1,1,0,0)
geography.matrix[[4]][,3]<-c(1,1,1,0,0)
geography.matrix[[4]][,4]<-c(0,0,0,1,1)
geography.matrix[[4]][,5]<-c(0,0,0,1,1)



cleanEx()
nameEx("InfTemp")
### * InfTemp

flush(stderr()); flush(stdout())

### Name: InfTemp
### Title: Paleotemperature data across the Cenozoic
### Aliases: InfTemp
### Keywords: datasets

### ** Examples

data(InfTemp)
plot(InfTemp)



cleanEx()
nameEx("JSDtree")
### * JSDtree

flush(stderr()); flush(stdout())

### Name: JSDtree
### Title: The Jensen-Shannon distance metric between phylogenies
### Aliases: JSDtree

### ** Examples

trees<-TESS::sim.globalBiDe.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
JSDtree(trees,alpha=0.8)



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
nameEx("Phyllostomidae")
### * Phyllostomidae

flush(stderr()); flush(stdout())

### Name: Phyllostomidae
### Title: Phyllostomidae phylogeny
### Aliases: Phyllostomidae
### Keywords: datasets

### ** Examples

data(Phyllostomidae)
print(Phyllostomidae)
plot(Phyllostomidae)



cleanEx()
nameEx("Phyllostomidae_genera")
### * Phyllostomidae_genera

flush(stderr()); flush(stdout())

### Name: Phyllostomidae_genera
### Title: Phylogenies of Phyllostomidae genera
### Aliases: Phyllostomidae_genera
### Keywords: datasets

### ** Examples

data(Phyllostomidae_genera)
print(Phyllostomidae_genera)



cleanEx()
nameEx("fitCompModel")
### * fitCompModel

flush(stderr()); flush(stdout())

### Name: fitCompModel
### Title: Fit models of trait evolution incorporating interspecific
###   interactions
### Aliases: fitCompModel

### ** Examples


data(Anolis.data)
geography.object<-Anolis.data$geography.object
pPC1<-Anolis.data$data
phylo<-Anolis.data$phylo

#Fit three models without biogeography to pPC1 data
MC.fit<-fitCompModel(phylo,pPC1,model="MC")
DDlin.fit<-fitCompModel(phylo,pPC1,model="DDlin")
DDexp.fit<-fitCompModel(phylo,pPC1,model="DDexp")

#Now fit models that incorporate biogeography, NOTE these models take longer to fit
MC.geo.fit<-fitCompModel(phylo,pPC1,model="MC",geography.object=geography.object)
DDlin.geo.fit<-fitCompModel(phylo,pPC1,model="DDlin",geography.object=geography.object)
DDexp.geo.fit<-fitCompModel(phylo,pPC1,model="DDexp",geography.object=geography.object)




cleanEx()
nameEx("fit_bd")
### * fit_bd

flush(stderr()); flush(stdout())

### Name: fit_bd
### Title: Maximum likelihood fit of the general birth-death model
### Aliases: fit_bd

### ** Examples

# Some examples may take a little bit of time. Be patient!

data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)

# Fit the pure birth model (no extinction) with a constant speciation rate
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.09)
mu_par<-c()
#result_cst <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                     f=87/89,cst.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
#result_cst$model <- "pure birth with constant speciation rate"

# Fit the pure birth model (no extinction) with exponential variation
# of the speciation rate with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
#result_exp <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                     f=87/89,expo.lamb=TRUE,fix.mu=TRUE,dt=1e-3)
#result_exp$model <- "pure birth with exponential variation in speciation rate"

# Fit the a pure birth model (no extinction) with linear variation of
# the speciation rate with time
f.lamb <-function(t,y){y[1] + y[2] * t}
f.mu<-function(t,y){0}
lamb_par<-c(0.09, 0.001)
mu_par<-c()
#result_lin <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89,fix.mu=TRUE,dt=1e-3)
#result_lin$model <- "pure birth with linear variation in speciation rate"

# Fit a birth-death model with exponential variation of the speciation
# rate with time and constant extinction
f.lamb<-function(t,y){y[1] * exp(y[2] * t)}
f.mu <-function(t,y){y[1]}
lamb_par <- c(0.05, 0.01)
mu_par <-c(0.005)
#result_bexp_dcst <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                           f=87/89,expo.lamb=TRUE,cst.mu=TRUE,dt=1e-3)
#result_bexp_dcst$model <- "birth-death with exponential variation in speciation rate
#                           and constant extinction"

# Find the best model
#index <- which.min(c(result_cst$aicc, result_exp$aicc, result_lin$aicc,result_bexp_dcst$aicc))
#rbind(result_exp, result_lin, result_cst, result_bexp_dcst)[index,]




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
nameEx("fit_env")
### * fit_env

flush(stderr()); flush(stdout())

### Name: fit_env
### Title: Maximum likelihood fit of the environmental birth-death model
### Aliases: fit_env

### ** Examples

data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

# Fits a model with lambda varying as an exponential function of temperature
# and mu fixed to 0 (no extinction).  Here t stands for time and x for temperature.
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
#result_exp <- fit_env(Cetacea,InfTemp,tot_time,f.lamb,f.mu,lamb_par,mu_par,
#                      f=87/89,fix.mu=TRUE,df=dof,dt=1e-3)



cleanEx()
nameEx("fit_sgd")
### * fit_sgd

flush(stderr()); flush(stdout())

### Name: fit_sgd
### Title: Maximum likelihood fit of the SGD model
### Aliases: fit_sgd

### ** Examples

# Some examples may take a little bit of time. Be patient!
data(Calomys)
tot_time <- max(node.age(Calomys)$ages)
par_init <- c(10, 1, 1)
#fit_sgd(Calomys, tot_time, par_init, f=11/13)



cleanEx()
nameEx("likelihoodDD")
### * likelihoodDD

flush(stderr()); flush(stdout())

### Name: likelihoodDD
### Title: Likelihood of a dataset under diversity-dependent models.
### Aliases: likelihoodDD

### ** Examples

data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data

# Compute the likelihood that the r value is twice the ML estimate for the DDexp model
par <- c(0.08148371, (2*-0.3223835))
lh <- -likelihoodDD(phylo,pPC1,par,model="DDexp")



cleanEx()
nameEx("likelihoodDD_geog")
### * likelihoodDD_geog

flush(stderr()); flush(stdout())

### Name: likelihoodDD_geog
### Title: Likelihood of a dataset under diversity-dependent models with
###   biogeography.
### Aliases: likelihoodDD_geog

### ** Examples

data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data
geography.object <- Anolis.data$geography.object

# Compute the likelihood with geography using ML parameters for fit without geography
par <- c(log(0.01153294),-0.0006692378)
lh <- -likelihoodDD_geog(phylo,pPC1,par,geography.object,model="DDlin")



cleanEx()
nameEx("likelihood_MC")
### * likelihood_MC

flush(stderr()); flush(stdout())

### Name: likelihood_MC
### Title: Likelihood of a dataset under the matching competition model.
### Aliases: likelihood_MC

### ** Examples

data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data

# Compute the likelihood that the S value is twice the ML estimate
par <- c(0.0003139751, (2*-0.06387258))
lh <- -likelihood_MC(phylo,pPC1,par)



cleanEx()
nameEx("likelihood_MC_geog")
### * likelihood_MC_geog

flush(stderr()); flush(stdout())

### Name: likelihood_MC_geog
### Title: Likelihood of a dataset under the matching competition model
###   with biogeography.
### Aliases: likelihood_MC_geog

### ** Examples

data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data
geography.object <-  Anolis.data$geography.object

# Compute the likelihood with geography using ML parameters for fit without geography
par <- c(0.0003139751, -0.06387258)
lh <- -likelihood_MC_geog(phylo,pPC1,par,geography.object)



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
lh <- likelihood_bd(Cetacea,tot_time,f.lamb,f.mu,f,cst.mu=TRUE,expo.lamb=TRUE, dt=1e-3)



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
likelihood <- likelihood_coal_cst(Vtimes,ntips,tau0,gamma,N0)



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



cleanEx()
nameEx("likelihood_sgd")
### * likelihood_sgd

flush(stderr()); flush(stdout())

### Name: likelihood_sgd
### Title: Likelihood of a phylogeny under the SGD model
### Aliases: likelihood_sgd

### ** Examples

data(Cetacea)
tot_time <- max(node.age(Cetacea)$ages)
b <- 1e6
d <- 1e6-0.5
nu <- 0.6
f <- 87/89
#lh <- likelihood_sgd(Cetacea, tot_time, b, d, nu, f)



cleanEx()
nameEx("plot_BICompare")
### * plot_BICompare

flush(stderr()); flush(stdout())

### Name: plot_BICompare
### Title: Display modalities on a phylogeny.
### Aliases: plot_BICompare

### ** Examples


data(Cetacea)
result <- BICompare(Cetacea,5)
#plot_BICompare(Cetacea,result)



cleanEx()
nameEx("plot_JSDtree")
### * plot_JSDtree

flush(stderr()); flush(stdout())

### Name: plot_JSDtree
### Title: Heatmap comparison of phylogenies
### Aliases: plot_JSDtree

### ** Examples

trees<-TESS::sim.globalBiDe.age(n=20,age=10,0.15,0.05,MRCA=TRUE)
res<-JSDtree(trees,alpha=0.8)
#plot_JSDtree(res)



cleanEx()
nameEx("plot_dtt")
### * plot_dtt

flush(stderr()); flush(stdout())

### Name: plot_dtt
### Title: Plot diversity through time
### Aliases: plot_dtt

### ** Examples


data(Balaenopteridae)
tot_time<-max(node.age(Balaenopteridae)$ages)

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate
# with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.08, 0.01)
mu_par<-c()
result <- fit_bd(Balaenopteridae,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=1,
                     expo.lamb = TRUE, fix.mu=TRUE)

# plot estimated number of species through time
# plot_dtt(result, tot_time, N0=9)



cleanEx()
nameEx("plot_fit_bd")
### * plot_fit_bd

flush(stderr()); flush(stdout())

### Name: plot_fit_bd
### Title: Plot speciation, extinction & net diversification rate functions
###   of a fitted model
### Aliases: plot_fit_bd

### ** Examples


data(Balaenopteridae)
tot_time<-max(node.age(Balaenopteridae)$ages)

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate
# with time
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.08, 0.01)
mu_par<-c()
result <- fit_bd(Balaenopteridae,tot_time,f.lamb,f.mu,lamb_par,mu_par,
                     expo.lamb = TRUE, fix.mu=TRUE)
# plot fitted rates
#plot_fit_bd(result, tot_time)



cleanEx()
nameEx("plot_fit_env")
### * plot_fit_env

flush(stderr()); flush(stdout())

### Name: plot_fit_env
### Title: Plot speciation, extinction & net diversification rate functions
###   of a fitted environmental model
### Aliases: plot_fit_env

### ** Examples


data(Balaenopteridae)
tot_time<-max(node.age(Balaenopteridae)$ages)
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

# Fit the pure birth model (no extinction) with exponential variation of the speciation rate
# with time
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()
result <- fit_env(Balaenopteridae,InfTemp,tot_time,f.lamb,f.mu,
      lamb_par,mu_par,f=1, fix.mu=TRUE, df=dof, dt=1e-3)

# plot fitted rates
#plot_fit_env(result, InfTemp, tot_time)



cleanEx()
nameEx("plot_spectR")
### * plot_spectR

flush(stderr()); flush(stdout())

### Name: plot_spectR
### Title: Spectral density plot of a phylogeny.
### Aliases: plot_spectR

### ** Examples


data(Cetacea)
result <- spectR(Cetacea)
#plot_spectR(result)



cleanEx()
nameEx("simCompModel")
### * simCompModel

flush(stderr()); flush(stdout())

### Name: simCompModel
### Title: Recursive simulation (root-to-tip) of competition models
### Aliases: simCompModel

### ** Examples


data(Cetacea)

# Simulate data under the matching competition model
MC.data<-simCompModel(Cetacea,pars=list(sig2=0.01,S=-0.1),root.value=0,Nsegments=1000,model="MC")

# Simulate data under the diversity dependent linear model
DDlin.data<-simCompModel(Cetacea,pars=list(sig2=0.01,b=-0.0001),root.value=0,Nsegments=1000,model="DDlin")

# Simulate data under the diversity dependent linear model
DDexp.data<-simCompModel(Cetacea,pars=list(sig2=0.01,r=-0.01),root.value=0,Nsegments=1000,model="DDexp")





cleanEx()
nameEx("sim_sgd")
### * sim_sgd

flush(stderr()); flush(stdout())

### Name: sim_sgd
### Title: Algorithm for simulating a phylogenetic tree under the SGD model
### Aliases: sim_sgd

### ** Examples

tau <- 10
b <- 1e6
d <- b-0.5
nu <- 0.6
tree <- sim_sgd(tau,b,d,nu)
plot(tree)



cleanEx()
nameEx("spectR")
### * spectR

flush(stderr()); flush(stdout())

### Name: spectR
### Title: Spectral density plot of a phylogeny
### Aliases: spectR

### ** Examples

data(Cetacea)
spectR(Cetacea,method="standard")



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
