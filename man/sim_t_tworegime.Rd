\name{sim_t_tworegime}
\alias{sim_t_tworegime}

\title{
Recursive simulation (root-to-tip) of two-regime models 
}
\description{
Simulates datasets for a given phylogeny under two-regime matching competition (MC), diversity dependent linear (DDlin), diversity dependent exponential (DDexp), or early burst (EB) models of trait evolution. Simulations are carried out from the root to the tip of the tree.
}
\usage{

sim_t_tworegime(regime.map, pars, root.value, Nsegments=2500, 
                model=c("MC","DDexp","DDlin","EB"),
	            	verbose=TRUE, rnd=6)

}

\arguments{
  \item{regime.map}{
a stochastic map of the two regimes stored as a simmap object output from \code{\link[phytools]{make.simmap}}
}
  \item{pars}{
a vector containing the three parameters for the chosen model; all models require \code{sig2}, and additionally, the MC model requires \code{S1} and \code{S2},
 specifying the level of competition in regime 1 and 2, respectively (larger negative values correspond to higher levels of competition), the DDlin model 
 requires \code{b1} and \code{b2}, the DDexp model requires \code{r1}, the slope parameters (negative in cases of decline in evolutionary 
 rates with increasing diversity). \code{sig2} must be listed first.
}
  \item{root.value}{
a number specifying the trait value for the ancestor
}
  \item{Nsegments}{
a value specifying the total number of time segments to simulate across for the phylogeny (see Details)
}
  \item{model}{
model chosen to fit trait data, \code{"MC"} is the matching competition model, \code{"DDlin"} is the diversity-dependent linear model,  \code{"DDexp"} is the diversity-dependent exponential model, and \code{"EB"} is the early burst model.
  }
  \item{verbose}{
if \code{TRUE}, prints the identity of regimes corresponding to each parameter value
}

  \item{rnd}{
  number of digits to round timings to (see \code{\link{round}} (see Details)
}


}

\value{
	a named vector with simulated trait values for \eqn{n} species in the phylogeny
}

\details{
Adjusting \code{Nsegments} will impact the length of time the simulations take. 
The length of each segment (\code{max(nodeHeights(phylo))/Nsegments}) 
should be much smaller than the smallest branch (\code{min(phylo$edge.length)}).

Adjusting \code{rnd} may help if function crashes.
}

\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
	
	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{J Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}}

\seealso{
\code{\link{fit_t_comp}}
}

\examples{

data(Cetacea_clades)

\donttest{

# Simulate data under the matching competition model
MC_tworegime.data<-sim_t_tworegime(Cetacea_clades,pars=c(sig2=0.01,S1=-0.1,S2=-0.01),
	root.value=0,Nsegments=1000,model="MC")

# Simulate data under the diversity dependent linear model
DDlin_tworegime.data<-sim_t_tworegime(Cetacea_clades,pars=c(sig2=0.01,b1=-0.0001,b2=-0.000001),
	root.value=0,Nsegments=1000,model="DDlin")

# Simulate data under the diversity dependent linear model
DDexp_tworegime.data<-sim_t_tworegime(Cetacea_clades,pars=c(sig2=0.01,r1=-0.01,r2=-0.02),
	root.value=0,Nsegments=1000,model="DDexp")

# Simulate data under the diversity dependent linear model
EB.data_tworegime<-sim_t_tworegime(Cetacea_clades,pars=c(sig2=0.01,r1=-0.01,r2=-0.02),
	root.value=0,Nsegments=1000,model="EB")


}

}
