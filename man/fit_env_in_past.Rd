\name{fit_env_in_past}
\alias{fit_env_in_past}

\title{
Maximum likelihood fit of the environmental birth-death model excluding the recent past
}
\description{
Fits the environmental birth-death model with potentially missing extant species to a phylogeny, by maximum likelihood while excluding the recent past. Notations follow Morlon et al. PNAS 2011 and Condamine et al. ELE 2013.
}
\usage{
fit_env_in_past(phylo, env_data, tot_time, time_stop, f.lamb, f.mu, desc, tot_desc, 
lamb_par, mu_par, df= NULL, meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE,
       expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
       dt=0, cond = "crown")
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation) that does not include any recent speciation (i.e. no speciation events between time_stop and the present). 
}
  \item{env_data}{
environmental data, given as a data frame with two columns. The first column is time, the second column is the environmental data (temperature for instance).
}
  \item{time_stop}{
the age of the phylogeny where to stop the birth-death process: it excludes the recent past (between the present and time_stop), while conditioning on the survival of the lineages from time_stop to the present.
}
  \item{tot_time}{
the age of the phylogeny (crown age, or stem age if known). If working with crown ages, tot_time is given by max(node.age(phylo)$ages).
}
  \item{f.lamb}{
a function specifying the hypothesized functional form of the variation of the speciation rate \eqn{\lambda} with time and the environmental variable. Any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the third arguement is a numeric vector of the parameters controlling the time and environmental variation (to be estimated).
}
  \item{f.mu}{
a function specifying the hypothesized functional form of the variation of the extinction rate \eqn{\mu} with time and the environmental variable. Any functional form may be used.
This function has three arguments: the first argument is time; the second argument is the environmental variable; the second argument is a numeric vector of the parameters controlling the time and environmental variation (to be estimated).
}
  \item{lamb_par}{
a numeric vector of initial values for the parameters of f.lamb to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong.
}
  \item{mu_par}{
a numeric vector of initial values for the parameters of f.mu to be estimated (these values are used by the optimization algorithm). The length of this vector is used to compute the total number of parameters in the model, so to fit a model without extinction (for example), mu_par should be empty (vector of length 0). Otherwise aic values will be wrong.
}
\item{df}{
 the degree of freedom to use to define the spline. As a default, smooth.spline(env_data[,1], env_data[,2])$df is used. See \emph{sm.spline} for details.  
}
  \item{desc}{
  the number of  lineages present at time_stop (i.e. the number of lineages present in the phylogenetic tree).

}
  \item{tot_desc}{
  the total number of extant species (included in the unsampled ones). 
}
  \item{meth}{
optimization to use to maximize the likelihood function, see \emph{optim} for more details.
}
  \item{cst.lamb}{
logical: should be set to TRUE only if f.lamb is constant (i.e. does not depend on time or the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{cst.mu}{
logical: should be set to TRUE only if f.mu is constant (i.e. does not depend on time or the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.lamb}{
logical: should be set to TRUE only if f.lamb is an exponential function of time (and does not depend on the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{expo.mu}{
logical: should be set to TRUE only if f.mu is an exponential function of time (and does not depend on the environmental variable) to use analytical instead of numerical computation in order to reduce computation time.
}
  \item{fix.mu}{
 logical: if set to TRUE, the extinction rate \eqn{\mu} is fixed and will not be optimized.
}
  \item{dt}{
the default value is 0. In this case, integrals in the likelihood are computed using R "integrate" function, which can be quite slow. If a positive dt is given as argument, integrals are computed using a piece-wise contant approximation, and dt represents the length of the intervals on which functions are assumed to be constant. We found that 1e-3 generally provides a good trade-off between precision and computation time.
}
  \item{cond}{
conditioning to use to fit the model:
\itemize{ \item FALSE: no conditioning (not recommended);
\item "stem": conditioning on the survival of the stem lineage (use when the stem age is known, in this case tot_time should be the stem age);
\item "crown" (default): conditioning on a speciation event at the crown age and survival of the 2 daugther lineages (use when the stem age is not known, in this case tot_time should be the crown age).
}
}
}


\value{
	a list with the following components
  \item{model}{the name of the fitted model}
  \item{LH}{the maximum log-likelihood value}
  \item{aicc}{the second order Akaike's Information Criterion}
  \item{lamb_par}{a numeric vector of estimated f.lamb parameters, in the same order as defined in f.lamb}
  \item{mu_par}{a numeric vector of estimated f.mu parameters, in the same order as defined in f.mu (if fix.mu is FALSE)}
}

\details{The lengths of lamb_par and mu_par are used to compute the total number of parameters in the model, so to fit a model with constant speciation rate (for example), lamb_par should be a vector of length 1. Otherwise aic values will be wrong. In the f.lamb and f.mu functions, time runs from the present to the past.
%The rate models here are function of environmental data fixed by
%env_data.
%Hence, if the parameter controlling the variation of \eqn{\lambda} with time is estimated to be positive (for example), this means that the speciation rate decreases from past to present.
Note that abs(f.lamb) and abs(f.mu) are used in the likelihood computation as speciation and extinction rates should always be positive. A consequence of this is that negative speciation/extinction rates estimates can be returned. They should be interpreted in aboslute terms. See Morlon et al. 2020 for a more detailed explanation.}

\note{
	The speed of convergence of the fit might depend on the degree of freedom chosen to define the spline. 
}
\references{
Morlon, H., Parsons, T.L. and Plotkin, J.B. (2011) Reconciling molecular phylogenies with the fossil record \emph{Proc Nat Acad Sci} 108: 16327-16332
	
Condamine, F.L., Rolland, J., and Morlon, H. (2013) Macroevolutionary perspectives to environmental change, Eco Lett 16: 72-85

	Lewitus, E., Bittner, L., Malviya, S., Bowler, C., & Morlon, H. (2018) Clade-specific diversification dynamics of marine diatoms since the Jurassic \emph{Nature Ecology and Evolution}, 2(11), 1715–1723
	
  Perez-Lamarque, B., Öpik, M.,Maliet, O., Afonso Silva, A., Selosse, M-A., Martos, F., Morlon, H., Global drivers of obligate mycorrhizal symbionts diversification, \emph{bioRxiv} doi: 10.1101/2020.07.28.224790
  
}

\author{H Morlon, F Condamine, E Lewitus, B Perez-Lamarque}

\seealso{
\code{\link{plot_fit_env}}, \code{\link{fit_bd_in_past}}, \code{\link{fit_env}}
}

\examples{

library(ape)
library(phytools)
library(pspline)

data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
data(InfTemp)
dof<-smooth.spline(InfTemp[,1], InfTemp[,2])$df

plot(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)

# slice the Cetaceae tree 10 Myr ago:
time_stop=10
sliced_tree <- Cetacea
sliced_sub_trees <- treeSlice(sliced_tree,slice = tot_time-time_stop, trivial=TRUE)
for (i in 1:length(sliced_sub_trees)){if (Ntip(sliced_sub_trees[[i]])>1){
  sliced_tree <- drop.tip(sliced_tree,tip=sliced_sub_trees[[i]]$tip.label[2:Ntip(sliced_sub_trees[[i]])])
}}
for (i in which(node.depth.edgelength(sliced_tree)>(tot_time-time_stop))){sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)] <- sliced_tree$edge.length[which(sliced_tree$edge[,2]==i)]-time_stop}

Ntip(sliced_tree) # 27 lineages present 10 Myr have survived until today

# Now we can fit environment-dependent birth-death models excluding the 10 last Myr

# Fits a model with lambda varying as an exponential function of temperature
# and mu fixed to 0 (no extinction).  Here t stands for time and x for temperature.
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.10, 0.01)
mu_par<-c()

#result_env <- fit_env_in_past(sliced_tree, InfTemp, tot_time, time_stop, f.lamb, 
#                             f.mu, lamb_par,mu_par,
#                             desc=Ntip(sliced_tree), tot_desc=89, 
#                             fix.mu=TRUE,df=dof,dt=1e-3)

}
