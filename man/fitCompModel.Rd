\name{simCompModel}
\alias{simCompModel}

\title{
Recursive simulation (root-to-tip) of competition models 
}
\description{
Simulates datasets for a given phylogeny under matching competition (MC), diversity dependent linear (DDlin), or diversity dependent exponential (DDexp) models of trait evolution. Simulations are carried out from the root to the tip of the tree.
}
\usage{

fitCompModel(phylo,data,model="MC,DDexp,DDlin",geography.object=NULL,method="Nelder-Mead")

}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{data}{
a named vector of trait values with names matching \code{phylo$tip.label}
}
  \item{model}{
model chosen to fit trait data, \code{"MC"} is the matching competition model of Nuismer & Harmon 2014, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}
  \item{geography.object}{
a list specifying sympatry (= 1) and allopatry (= 0) for each extant lineages
at each internode interval (see details)
}
  \item{method}{
method to use for ML optimization (default is Nelder-Mead)
}


\value{
	a list with the following elements:
	
}

\details{
Adjusting \code{\Nsegments} will impact the length of time the simulations take. 
The length of each segment (max(nodeHeights(phylo))/Nsegments) 
should be much smaller than the smallest branch (min(phylo$edge.length)).

\code{\geography.object} should be a list of matrices with row and column
names that match the "nat" object created within the VCV.rescale() function.

}

\references{
	Drury, J., Clavel, J., and Morlon, H. (submitted) Estimating the effect of competition on trait evolution using maximum likelihood inference.

	Nuismer & Harmon
	
	Weir & Mursleen
}

\author{Jonathan Drury \href{mailto:drury@biologie.ens.fr}{drury@biologie.ens.fr}
		Julien Clavel
}

\seealso{
\code{\link{simCompModel}}
\code{\link{CreateGeoObject}}
}

\note{
In current version, S parameter is restricted to take on negative values in MC + geography ML optimization.
}
\examples{

data(Calomys)
data(Calomys.geo)

# First, simulate data under the matching competition model
MC.data<-simCompModel(Calomys,pars=list(sig2=0.05,S=-1),root.value=0,Nsegments=1000,model="MC")

#Now fit three models to these data
MC.fit<-fitCompModel(Calomys,MC.data,model="MC")
DDlin.fit<-fitCompModel(Calomys,MC.data,model="DDlin")
DDexp.fit<-fitCompModel(Calomys,MC.data,model="DDlin")

#Now fit models that incorporate biogeography (to see form of geography.object, examine Calomys.geo)
##NOTE: Calomys.geo is purely illustrative and not based on actual data
MC.geo.fit<-fitCompModel(Calomys,MC.data,model="MC",geography.object=Calomys.geo)
DDlin.geo.fit<-fitCompModel(Calomys,MC.data,model="DDlin",geography.object=Calomys.geo)
DDexp.geo.fit<-fitCompModel(Calomys,MC.data,model="DDexp",geography.object=Calomys.geo)


}
