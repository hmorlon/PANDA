\name{Anolis.data}
\alias{Anolis.data}

\title{
Anolis dataset
}
\description{
Phylogeny, data, and geography.object for a subclade of Greater Antillean \emph{Anolis} lizards. 

}
\usage{
data(Anolis.data)
}

\details{
	Phylogeny trimmed from the maximum clade credibility tree of Mahler et al. 2013, phylogenetic principal component data are from Mahler et al. 2013, and biogeography data are from Mahler & Ingram 2014. 
}

\references{
	Drury, J., Clavel, J., and Morlon, H. (submitted) Estimating the effect of competition on trait evolution using maximum likelihood inference.

	Mahler, D.L., Ingram, T., Revell, L., and Losos, J. 2013. Exceptional convergence on the macroevolutionary landscape in island lizard radiations. \emph{Science}. 341:292-295.
	
	Mahler, D.L. and Ingram, T. 2014. Phylogenetic comparative methods for studying clade-wide convergence. In \emph{Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology}, ed. L. Garamszegi. pp.425-450.
	}


\examples{
data(Anolis.data)
plot(Anolis.data$phylo)
print(Anolis.data$data)
print(Anolis.data$geography.object[[16]])

# Compute the likelihood that the S value is twice the ML estimate
par <- c(0.0003139751, (2*-0.06387258))
lh <- -likelihood_MC(phylo,pPC1,par)
}
