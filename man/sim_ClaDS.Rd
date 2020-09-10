\name{sim_ClaDS}
\alias{sim_ClaDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the ClaDS model
}
\description{
Simulate a birth-death phyloh-geny with rate shifts happening at speciation events.
}
\usage{
sim_ClaDS(lambda_0, mu_0,
          new_lamb_law="lognormal*shift",new_mu_law="turnover",
          condition="time", time_stop = 0, taxa_stop = Inf,
          sigma_lamb=0.1, alpha_lamb=1, lamb_max=1,lamb_min=0,
          sigma_mu=0, alpha_mu=1, mu_min=mu_0,mu_max=mu_0, 
          theta=1,nShiftMax=Inf,
          return_all_extinct=FALSE,prune_extinct=TRUE,
          maxRate=Inf)
                       }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda_0}{
Initial speciation rate.
}
  \item{mu_0}{
Initial extinction rate, or turnover rate if new_mu_law == "turnover".
}
  \item{new_lamb_law}{
Distribution in which the new speciation rates are drawn at a speciation event. See details.
}
  \item{new_mu_law}{
Distribution in which the new extinction rates are drawn at a speciation event. See details.
}
  \item{condition}{
Stoping condition. Can be "time" (the default) or "taxa".
}
  \item{time_stop}{
Stoping time if condition == "time".
}
  \item{taxa_stop}{
Final number of species if condition == "taxa". 

If condition == "time", the process is stoped if the number of species exceeds taxa_stop. This can be usefull for some parametrizations of the model for which the number of species can reach very large number very quickly, leading to computation time and memory issues. To disable this option, use taxa_stop = Inf (the default).
}
  \item{sigma_lamb}{
Parameter of the new speciation rates distribution, see details.
}
  \item{alpha_lamb}{
Parameter of the new speciation rates distribution, see details.
}
  \item{lamb_max}{
Parameter of the new speciation rates distribution, see details.
}
  \item{lamb_min}{
Parameter of the new speciation rates distribution, see details.
}
  \item{sigma_mu}{
Parameter of the new extinction rates distribution, see details.
}
  \item{alpha_mu}{
Parameter of the new extinction rates distribution, see details.
}
  \item{mu_min}{
Parameter of the new extinction rates distribution, see details.
}
  \item{mu_max}{
Parameter of the new extinction rates distribution, see details.
}
\item{theta}{
Probability to have a rate shift at speciation. Default to 1.
}
  \item{nShiftMax}{
Maximum number of rate shifts. If nShiftMax < Inf, theta is set to 0 as soon as there has been nShiftMax rate shifts. Set nShiftMax = Inf (the default) to disable this option.
}
  \item{return_all_extinct}{
Boolean specifying whether the function should return extinct phylogenies. Default to FALSE.
}
  \item{prune_extinct}{
Boolean specifying whether extinct species should be removed from the resulting phylogeny. Default to TRUE.
}
  \item{maxRate}{
The process is stoped if one of the lineage has a speciation rate that exceeds maxRate. This can be usefull for some parametrizations of the model for which the rates can reach very large values, leading to numerical overflows. To disable this option, use maxRate = Inf (the default).
}
}
\details{
Available options for new_lamb_law are : 

\itemize{
\item "uniform", the new speciation rates are drawn uniformly in [lamb_min, lamb_max].
\item "normal", the new speciation rates are drawn in a normal distribution with parameters (sigma_lamb^2, parent_lambda), truncated in 0.
\item "lognormal", the new speciation rates are drawn in a lognormal distribution with parameters (sigma_lamb^2, parent_lambda).
\item "lognormal*shift", the new speciation rates are drawn in a lognormal distribution with parameters (sigma_lamb^2, parent_lambda * alpha_lamb). This is the default option as it corresponds to the ClaDS model.
\item "lognormal*t", the new speciation rates are drawn in a lognormal distribution with parameters (sigma_lamb^2 * t^2, parent_lambda), where t is the age of the mother species.
\item "logbrownian", the new speciation rates are drawn in a lognormal distribution with parameters (sigma_lamb^2 * t, parent_lambda), where t is the age of the mother species. This is used to approximate the case where speciation rates are evolving as the log of a brownian motion, as is done in Beaulieu, J. M. and B. C. O'Meara. (2015).
\item "normal+shift", the new speciation rates are drawn in a normal distribution with parameters (sigma_lamb^2, parent_lambda + alpha_lamb), truncated in 0.
\item "normal*shift", the new speciation rates are drawn in a normal distribution with parameters (sigma_lamb^2, parent_lambda * alpha_lamb), truncated in 0.
}

Available options for new_mu_law are : 

\itemize{
\item "uniform", the new extinction rates are drawn uniformly in [mu_min, mu_max].
\item "normal", the new extinction rates are drawn in a normal distribution with parameters (sigma_mu^2, parent_mu), truncated in 0.
\item "lognormal", the new extinction rates are drawn in a lognormal distribution with parameters (sigma_mu^2, parent_mu).
\item "lognormal*shift", the new extinction rates are drawn in a lognormal distribution with parameters (sigma_mu^2, parent_mu * alpha_mu).
\item "normal*t", the new speciation rates are drawn in a normal distribution with parameters (sigma_lamb^2 * t^2, parent_lambda), where t is the age of the mother species.
\item "turnover", the turnover rate is constant (in that case mu_0 is the turnover rate), so the new extinction rates are mu_0 times the new speciation rates. This is the default option, corresponding to ClaDS2.
}
}
\value{
A list with :
\item{tree }{The resulting phylogeny.}
\item{times }{A vector with the times of all speciation and extinction events.}
\item{nblineages }{A vector in which nblineages[i] is the number of species in the clade after the event happening at time times[i].}
\item{lamb }{A vector with all the different speciation rates resulting from the simulation.}
\item{mu }{A vector with all the different extinction rates resulting from the simulation.}
\item{rates }{A vector of integer mapping the elements of .$lamb and .$mu to the branches of .$tree.}
\item{maxRate }{A boolean indicating whether the process was ended before reaching the specified stopping criterion because one of the speciation rates exceeded maxRate (see the "arguments" section).}
\item{root_length }{The time before the first speciation event.}

}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution}, doi 10.1038/s41559-019-0908-0

Beaulieu, J. M. and B. C. O'Meara. 2015. Extinction can be estimated from moderately sized
molecular phylogenies. \emph{Evolution} 69:1036-1043.
}

\author{
O. Maliet
}
\seealso{
\code{\link{plot_ClaDS_phylo}}
}
\examples{
# Simulation of a ClaDS2 phylogeny
set.seed(1)

obj= sim_ClaDS( lambda_0=0.1,    
                mu_0=0.5,      
                sigma_lamb=0.7,         
                alpha_lamb=0.90,     
                condition="taxa",    
                taxa_stop = 20,    
                prune_extinct = TRUE)  

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

plot_ClaDS_phylo(tree,speciation_rates)


# Simulation of a phylogeny with constant extinction rate and speciation 
# rates evolving as a logbrownian
set.seed(4321)

obj= sim_ClaDS( lambda_0=0.1,    
                mu_0=0.2,    
                new_mu_law = "uniform",
                new_lamb_law = "logbrownian",
                sigma_lamb=0.4,         
                condition="taxa",    
                taxa_stop = 20,    
                prune_extinct = FALSE)  

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

par(mar=c(1,1,0,0))
plot_ClaDS_phylo(tree,speciation_rates)



# Simulation of a phylogeny with constant extinction rate and at most one shift
# in speciation rates
set.seed(1221)

obj= sim_ClaDS( lambda_0=0.1,    
                mu_0=0.05,    
                new_mu_law = "uniform",
                new_lamb_law = "uniform",
                lamb_max = 0.5, lamb_min = 0,     
                theta = 0.1, nShiftMax = 1,
                condition="taxa",    
                taxa_stop = 100,    
                prune_extinct = TRUE)  

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

plot_ClaDS_phylo(tree,speciation_rates)
}
