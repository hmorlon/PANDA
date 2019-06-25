\name{Caprimulgidae_ClaDS2}
\alias{Caprimulgidae_ClaDS2}
\docType{data}
\title{
An example run of ClaDS2.
}
\description{
An example of the run on the inference of ClaDS2 on the \emph{Caprimulgidae} phylogeny, thinned every 10 iterations.
}
\usage{data("Caprimulgidae_ClaDS2")}
\format{
A list object with fields :

\describe{
\item{\code{tree}}{The \emph{Caprimulgidae} phylogeny on which we ran the model.}
\item{\code{sample_fraction}}{The sample fraction for the clade.}
\item{\code{sampler}}{The chains obtained by running ClaDS2 on the \emph{Caprimulgidae} phylogeny.}
}
}
\details{
The \emph{Caprimulgidae} phylogeny was obtained from Jetz et al. (2012)
}
\source{
Jetz, W., G. Thomas, J. Joy, K. Hartmann, and A. Mooers. 2012. The global diversity of birds in
space and time. \emph{Nature 491:444.}
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution}, doi 10.1038/s41559-019-0908-0
}
\author{
O. Maliet
}
\seealso{
  \code{\link{fit_ClaDS}}, \code{\link{plot_ClaDS_chains}}, \code{\link{getMAPS_ClaDS0}}
}
\examples{
data("Caprimulgidae_ClaDS2")

# plot the mcmc chains
plot_ClaDS_chains(Caprimulgidae_ClaDS2$sampler)

\donttest{
# extract the Maxima A Posteriori for each parameter
maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1)
print(paste0("sigma = ", maps[1], " ; alpha = ", 
  maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))
  
# plot the infered branch specific speciation rates
plot_ClaDS_phylo(Caprimulgidae_ClaDS2$tree, maps[-(1:4)])
}
}