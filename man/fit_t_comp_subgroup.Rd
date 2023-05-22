\name{fit_t_comp_subgroup}
\alias{fit_t_comp_subgroup}

\title{
Fits models of trait evolution incorporating competitive interactions, restricting competition to occur only between members of a subgroup
}
\description{
Fits matching competition (MC), diversity dependent linear (DDlin), or diversity dependent exponential (DDexp) models of trait evolution to a given dataset, phylogeny, and stochastic maps of both subgroup membership and biogeography.
}
\usage{

fit_t_comp_subgroup(full.phylo, data, subgroup, subgroup.map, model=c("MC","DDexp","DDlin"),
	ana.events=NULL, clado.events=NULL, stratified=FALSE, regime.map=NULL,error=NULL, 
	par=NULL, method="Nelder-Mead", bounds=NULL)
	
}

\arguments{
  \item{full.phylo}{
an object of type 'phylo' (see ape documentation) containing all of the tips used to estimate ancestral biogeography in BioGeoBEARS
}
  \item{data}{
a named vector of trait values for subgroup members with names matching \code{full.phylo$tip.label}
}
  \item{subgroup}{
  subgroup whose members are competing
}
  \item{subgroup.map}{
a phylo object created using \code{\link{make.simmap}} in phytools that contains reconstructed subgroup membership
}
  \item{model}{
model chosen to fit trait data, \code{"MC"} is the matching competition model of Nuismer & Harmon 2014, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}

  \item{ana.events}{
the "ana.events" table produced in BioGeoBEARS that lists anagenetic events in the stochastic map
}
  \item{clado.events}{
the "clado.events" table produced in BioGeoBEARS that lists cladogenetic events in the stochastic map
}
  \item{stratified}{
logical indicating whether the stochastic map was built from a stratified analysis in BioGeoBEARS
}
 \item{error}{
A named vector with standard error (SE) for each species (with names matching \code{"phylo$tip.label"}). Default is NULL, if NA, then the SE is estimated from the data (a nuisance parameter for unknown errors). Note: When standard error are provided the nuisance parameter is also estimated.
}
  \item{par}{
vector specifying starting parameter values for maximum likelihood optimization. If unspecified, default values are used (see Details)
}
  \item{method}{
optimization algorithm to use (see \code{\link{optim}}; for DD models without biogeography, \code{\link{BB}} is also supported)
}
  \item{bounds}{
(optional) list of bounds to pass to optimization algorithm (see details at \code{\link{optim}})
}

}

\value{
	a list with the following elements:
  \item{LH}{maximum log-likelihood value}
  \item{aic}{Akaike Information Criterion value}
  \item{aicc}{AIC value corrected for small sample size}
  \item{free.parameters}{number of free parameters from the model}
  \item{sig2}{maximum-likelihood estimate of \code{sig2} parameter}
  \item{S}{maximum-likelihood estimate of \code{S} parameter of matching competition model (see Note)}
  \item{b}{maximum-likelihood estimate of \code{b} parameter of linear diversity dependence model (see Note)}
  \item{r}{maximum-likelihood estimate of \code{r} parameter of exponential diversity dependence model (see Note)}
  \item{z0}{maximum-likelihood estimate of \code{z0}, the value at the root of the tree}
  \item{convergence}{convergence diagnostics from \code{\link{optim}} function (see optim documentation)}
  \item{nuisance}{maximum-likelihood estimate of \code{nuisance}, the unknown, nuisance contribution to measurement error when \code{error} argument is used (that is NA or a vector provided by the user)}
}

\details{
If unspecified, \code{par} takes the default values of \code{var(data)/max(nodeHeights(phylo))} for sig2 and 0 for either \code{S} for the matching competition model, 
\code{b} for the linear diversity dependence model, or \code{r} for the exponential diversity dependence model. Values can be manually entered as a vector with the first element
equal to the desired starting value for sig2 and the second value equal to the desired starting value for either \code{S}, \code{b}, or \code{r}. Note: since likelihood optimization
uses sig rather than sig2, and since the starting value for is exponentiated to stabilize the likelihood search, if you input a \code{par} value, the first value specifying sig2 
should be the log(sqrt()) of the desired sig2 starting value. We recommend running ML optimization with several different starting values to ensure convergence.

Currently, this function can be used to implement the following models:
1. Subgroup pruning with biogeography: matching competition, diversity dependent
2. Subgroup pruning without biogeography: diversity dependent
3. Subgroup pruning without biogeography (two-regimes): diversity dependent (for more details, see \code{\link{fit_t_comp}}

}

\references{
	
	Drury, J., Clavel, J. Tobias, J., Rolland, J., Sheard, C., and Morlon, H. Tempo and mode of morphological evolution are decoupled from latitude in birds. \emph{PLOS Biology} doi:10.1371/journal.pbio.3001270

	Drury, J., Tobias, J., Burns, K., Mason, N., Shultz, A., and Morlon, H. 2018. Contrasting impacts of competition on ecological and social trait evolution in songbirds. \emph{PLOS Biology} 16(1): e2003563.

	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} 65: 700-710

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
	
	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

}

\seealso{
\code{\link{likelihood_subgroup_model}}
\code{\link{CreateGeobyClassObject}}
\code{\link{fit_t_comp}}

}

\note{
In current version, the \code{S} parameter is restricted to take on negative values in MC + geography ML optimization.
}
\examples{


\donttest{

data(BGB.examples)

#Prepare dataset with subgroups and biogeography

Canidae.phylo<-BGB.examples$Canidae.phylo
dummy.group<-c(rep("B",3),rep("A",12),rep("B",2),rep("A",6),rep("B",5),rep("A",6))
names(dummy.group)<-Canidae.phylo$tip.label


Canidae.simmap<-phytools::make.simmap(Canidae.phylo,dummy.group)

set.seed(123)
Canidae.data<-rnorm(length(Canidae.phylo$tip.label))
names(Canidae.data)<-Canidae.phylo$tip.label
Canidae.A<-Canidae.data[which(dummy.group=="A")]


#Fit model with subgroup pruning and biogeography
MC.fit_subgroup_geo<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events, 
clado.events=BGB.examples$Canidae.clado.events,stratified=FALSE,subgroup.map=Canidae.simmap, 
data=Canidae.A,subgroup="A",model="MC")

DDexp.fit_subgroup_geo<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events, 
clado.events=BGB.examples$Canidae.clado.events,stratified=FALSE,subgroup.map=Canidae.simmap, 
data=Canidae.A,subgroup="A",model="DDexp")

DDlin.fit_subgroup_geo<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events, 
clado.events=BGB.examples$Canidae.clado.events,stratified=FALSE,subgroup.map=Canidae.simmap, 
data=Canidae.A,subgroup="A",model="DDlin")

#Fit model with subgroup pruning and no biogeography (for DD models only)
DDexp.fit_subgroup_no.geo<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,data=Canidae.A,subgroup="A",
subgroup.map=Canidae.simmap,model="DDexp")

DDlin.fit_subgroup_no.geo<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,data=Canidae.A,subgroup="A",
subgroup.map=Canidae.simmap,model="DDlin")


#Prepare regime map for fitting two-regime models with subgroup pruning (for DD models only)
regime<-c(rep("regime1",15),rep("regime2",19))
names(regime)<-Canidae.phylo$tip.label
regime.map<-phytools::make.simmap(Canidae.phylo,regime)

#Fit model with subgroup pruning and two-regimes (for DD models only)
DDexp.fit_subgroup_two.regime<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,data=Canidae.A,subgroup="A",
subgroup.map=Canidae.simmap,model="DDexp",regime.map=regime.map)

DDlin.fit_subgroup_two.regime<-fit_t_comp_subgroup(full.phylo=Canidae.phylo,data=Canidae.A,subgroup="A",
subgroup.map=Canidae.simmap,model="DDlin",regime.map=regime.map)

	
}

}
