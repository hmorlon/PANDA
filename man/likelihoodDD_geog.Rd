\name{likelihoodDD_geog}
\alias{likelihoodDD_geog}

\title{
Likelihood of a dataset under diversity-dependent models with biogeography.
}
\description{
Computes the likelihood of a dataset under either the linear or exponential diversity dependent model with specified \code{sigma2} and slope values and with a \code{geography.object} formed using \code{\link{CreateGeoObject}}.  
}

\usage{
likelihoodDD_geog(phylo, data, par=c(sig2,slope),geography.object,model=c("DDlin","DDexp"))
}

\arguments{
  \item{phylo}{  
an object of type 'phylo' (see ape documentation)
}
  
  \item{data}{
a named vector of continuous data with names corresponding to \code{phylo$tip.label}
}

 \item{par}{
a vector listing a value for \code{sig2} and either \code{b} (for the linear diversity dependent model) or \code{r} (for the exponential diversity dependent model).
}

  \item{geography.object}{
a list specifying sympatry (= 1) and allopatry (= 0) for each extant lineages at each internode interval (see details)}

 \item{model}{
model chosen to fit trait data, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.


}

}

\details{
	When specifying \code{par}, \code{sig2} must be listed before the slope parameter (\code{b} or \code{r}).

	The \code{geography.object} should be a list of matrices with row and column
names that match the \code{$geography.object} object created with the \code{\link{CreateGeoObject}} function.

}

\note{
	To stabilize optimization, this function exponentiates the input \code{sig2} value, thus the user must input the log(sig2) value to compute the correct log likelihood (see example).
}	

\value{
  the negative log-likelihood value of the dataset (accordingly, the negative of the output should be recorded as the likelihood), given the phylogeny,  \code{sig2} and slope values, and \code{geography.object}.
}

\references{
	Drury, J., Clavel, J., and Morlon, H. (submitted) Estimating the effect of competition on trait evolution using maximum likelihood inference.

	Weir, J. & Mursleen, S. 2012. Diversity-dependent cladogenesis and trait evolution in the adaptive radiation of the auks (Aves: Alcidae). \emph{Evolution} 67:403-416.
}

\author{Jonathan Drury \href{mailto:drury@biologie.ens.fr}{drury@biologie.ens.fr}

		Julien Clavel
}

\seealso{
\code{\link{fitCompModel}}
\code{\link{CreateGeoObject}}
\code{\link{likelihoodDD}}

}


\examples{
data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data
geography.object <- Anolis.data$geography.object

# Compute the likelihood with geography using ML parameters for fit without geography
par <- c(log(0.01153294),-0.0006692378)
lh <- -likelihoodDD_geog(phylo,pPC1,par,geography.object,model="DDlin")
}
