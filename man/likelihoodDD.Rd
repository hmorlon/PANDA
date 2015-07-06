\name{likelihoodDD}
\alias{likelihoodDD}

\title{
Likelihood of a dataset under diversity-dependent models.
}
\description{
Computes the likelihood of a dataset under either the linear or exponential diversity dependent model with specified \code{sigma2} and slope values. 
}

\usage{
likelihoodDD(phylo, data, par=c(sig2,slope),model=c("DDlin","DDexp"))
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

 \item{model}{
model chosen to fit trait data, \code{"DDlin"} is the diversity-dependent linear model, and \code{"DDexp"} is the diversity-dependent exponential model of Weir & Mursleen 2013.
}

}

\details{
	When specifying \code{par}, \code{sig2} must be listed before the slope parameter (\code{b} or \code{r}).
}

\value{
  the negative log-likelihood value of the dataset (accordingly, the negative of the output should be recorded as the likelihood), given the phylogeny and \code{sig2} and slope values
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
\code{\link{likelihoodDD_geog}}

}


\examples{
data(Anolis.data)
phylo <- Anolis.data$phylo
pPC1 <- Anolis.data$data

# Compute the likelihood that the r value is twice the ML estimate for the DDexp model
par <- c(0.08148371, (2*-0.3223835))
lh <- -likelihoodDD(phylo,pPC1,par,model="DDexp")
}
