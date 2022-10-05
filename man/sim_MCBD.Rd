\name{sim_MCBD}
\alias{sim_MCBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of macroevolutionary diversification under the integrated model described in Aristide & Morlon 2019
}
\description{
Simulates the joint diversification of species and a continuous trait, where changes in both dimensions are interlinked through competitive interactions.
}
\usage{
sim_MCBD(pars, root.value = 0, age.max = 50, step.size = 0.01, bounds = c(-Inf,Inf),
         plot = TRUE, ylims=NULL, full.sim = FALSE)
}

\arguments{
  \item{pars}{
    Vector of simulation parameters:
    
    \code{pars[1]} corresponds to \eqn{lambda1}, the speciation intitation rate
    
    \code{pars[2]} corresponds to \eqn{tau0}, the basal speciation completion rate
    
    \code{pars[3]} corresponds to \eqn{beta}, the effect of trait differences on the speciation completion rate
    
    \code{pars[4]} corresponds to \eqn{mu0}, the competitive extinction parameter for good species
    
    \code{pars[5]} corresponds to \eqn{mubg}, the background good species extinction rate
    
    \code{pars[6]} corresponds to \eqn{mui0}, the competitive extinction parameter for incipient species
    
    \code{pars[7]} corresponds to \eqn{muibg}, the background incipient species extinction rate
    
    \code{pars[8]} corresponds to \eqn{alpha1}, the competition effect on extinction (competition strength)
    
    \code{pars[9]} corresponds to \eqn{alpha1}, the competition effect on trait evolution (competition strength)
    
    \code{pars[10]} corresponds to \eqn{sig2}, the variance (rate) of the Brownian motion
    
    \code{pars[11]} corresponds to \eqn{m}, the relative contribution of character displacement (competition) with respect to stochastic (brownian) evolution
}
\item{root.value}{the starting trait value}
\item{age.max}{maximum time for the simulation (if the process doesn't go extinct)}
\item{step.size}{size of each simulation step}
\item{bounds}{lower and upper value for bounds in trait space}
\item{plot}{logical indicating wether to plot the simulation}
\item{ylims}{y axis (trait values) limits for the simulation plot}
\item{full.sim}{logical indicating wether to return the full simulation (see details)}

}
\details{
It might be difficult to find parameter combinations that are sensitive. It is recommended to use the parameter settings of the examples as a staring point and from there modify them to understand the behaviour of the model.
If trees produced are too big, simulation can become too slow to ever finish. 
}
\value{
returns a list with the following elements:

\bold{all} contains the complete tree of the process (extant and extinct good and incipient lineages) and trait values for each tip in the tree

\bold{gsp_fossil} contains the extant and extinct good species tree and trait values for each tip in the tree

\bold{gsp_extant} contains the reconstructed (extant only) good species tree and trait values for each tip in the tree


If full.sim = TRUE, two additional elements are returned inside \bold{all}:

note: both elements are used internally to keep track of the simulation and are dynamically updated, so returned elements only reflect the last state
 
\bold{lin_mat} a matrix with information about the diversification process. Each row represents a new lineage in the process with the following elements:
  - Parental node, descendent node (0 if a tip), starting time, ending time, status at end (extinct(-2); incipient(-1); good(1)), speciation completion or extinction time; speciation completion time (NA if still incipient).

\bold{trait_mat} a list with trait values for each lineage at each time step throghout the simulation. Each element is a vector composed of the following:
Lineage number (same as row number in lin_mat), status (as in lin_mat), sister lineage number, trait values (NA if lineage didn't exist yet at that time step)


}
\author{
Leandro Aristide (leandroaristi@gmail.com)
}
\references{
Aristide, L., and Morlon, H. 2019. Understanding the effect of competition during evolutionary radiations: an integrated model of phenotypic and species diversification
}
\examples{
lambda1 = 0.25
tau0 = 0.01
beta = 0.6
mu0 = 0.5
mubg = 0.01
mui0 = 0.8
muibg = 0.02
alpha1 = alpha2 = 0.04
sig2 = 0.5
m = 20

pars <- c(lambda1, tau0, beta, mu0, mubg,mui0, muibg, alpha1, alpha2, sig2, m)

\dontshow{test = FALSE}
if(test){

#1000 steps, unbounded
res <- sim_MCBD(pars, age.max=10, step.size=0.01) 

#asymmetric bounds
res <- sim_MCBD(pars, age.max=10, step.size=0.01, bounds=c(-10,Inf)) 

#only deterministic component
pars <- c(lambda1, tau0, beta, mu0, mubg, mui0, muibg, alpha1, alpha2, sig2=0, m)
res <- sim_MCBD(pars, age.max=10)

plot(res$gsp_extant$tree)

}

}

