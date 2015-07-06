\name{fitCompModel}
\alias{fitCompModel}

\title{
Fit models of trait evolution incorporating interspecific interactions
}
\description{
Fits matching competition (MC), diversity dependent linear (DDlin), or diversity dependent exponential (DDexp) models of trait evolution to a given dataset and phylogeny.
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
}

\value{
	a list with the following elements:
  \item{lnL}{log-likelihood of maximum-likelihood fit}
  \item{sig2}{maximum-likelihood estimate of \code{sig2} parameter}
  \item{S}{maximum-likelihood estimate of \code{S} parameter of matching competition model (see Note)}
  \item{b}{maximum-likelihood estimate of \code{b} parameter of linear diversity dependence model}
  \item{r}{maximum-likelihood estimate of \code{r} parameter of exponential diversity dependence model}
  \item{aic}{Akaike Information Criterion value}
  \item{aicc}{AIC value corrected for small sample size}
  \item{free.parameter}{number of free parameters from the model}
  \item{convergence}{convergence diagnostics from \code{\link{optim}} function (see optim documentation)}

}

\details{

\code{geography.object} should be a list of matrices with row and column
names that match the \code{$geography.object} object created within the \code{\link{CreateGeoObject}} function.

}

\references{
	Drury, J., Clavel, J., and Morlon, H. (submitted) Estimating the effect of competition on trait evolution using maximum likelihood inference.

	Nuismer, S. & Harmon, L. 2015. Predicting rates of interspecific interaction from phylogenetic trees. \emph{Ecology Letters} 18:17-27.
	
	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{Jonathan Drury \href{mailto:drury@biologie.ens.fr}{drury@biologie.ens.fr}

		Julien Clavel
}

\seealso{
\code{\link{simCompModel}}
\code{\link{CreateGeoObject}}
}

\note{
In current version, the \code{S} parameter is restricted to take on negative values in MC + geography ML optimization.
}
\examples{

data(Anolis.data)
geography.object<-Anolis.data$geography.object
pPC1<-Anolis.data$data
phylo<-Anolis.data$phylo

#Fit three models without biogeography to pPC1 data
\code{MC.fit<-fitCompModel(phylo,pPC1,model="MC")}
\code{DDlin.fit<-fitCompModel(phylo,pPC1,model="DDlin")}
\code{DDexp.fit<-fitCompModel(phylo,pPC1,model="DDexp")}

#Now fit models that incorporate biogeography, NOTE these models take longer to fit
\code{MC.geo.fit<-fitCompModel(phylo,pPC1,model="MC",geography.object=geography.object)}
\code{DDlin.geo.fit<-fitCompModel(phylo,pPC1,model="DDlin",geography.object=geography.object)}
\code{DDexp.geo.fit<-fitCompModel(phylo,pPC1,model="DDexp",geography.object=geography.object)}

}
