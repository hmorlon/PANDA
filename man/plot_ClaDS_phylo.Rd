\name{plot_ClaDS_phylo}
\alias{plot_ClaDS_phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a phylogeny with branch-specific values
}
\description{
Plot a phylogeny with branches colored according to branch-specific rate values
}
\usage{
plot_ClaDS_phylo(phylo, rates, rates2 = NULL, 
                same.scale = T, main = NULL, lwd = 2, log = T, show.tip.label = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylo}{An object of class 'phylo'.}
  \item{rates}{A vector containing the branch-specific rates, in the same order as phylo$edges.}
  \item{rates2}{An optional second vector containing the branch-specific rates, in the same order as phylo$edges. If NULL (the default), the tree is only plotted once with the rate values from rates. If not, the tree is plotted twice, with the rate values from rates in the left panel and those from rates2 in the right panel. }
  \item{same.scale}{A boolean specifying whether the values from rates and rates2 are plotted with the same colorscale. Default to TRUE.}
  \item{main}{A title for the plot.}
  \item{lwd}{Width of the tree branch lengths. Default to 2.}
  \item{log}{A boolean specifying whether the rates values are plotted on a log scale. Default to TRUE.}
  \item{show.tip.label}{A boolean specifying whether the labels of the phylogeny should be displayed. Default to FALSE.}
  \item{...}{Optional arguments for \code{plot.phylo}.}
}
\references{
Maliet O., Hartig F. and Morlon H. 2019, A model with many small shifts for estimating species-specific diversificaton rates, \emph{Nature Ecology and Evolution}, doi 10.1038/s41559-019-0908-0
}
\author{
O. Maliet}
\seealso{
\code{\link{sim_ClaDS}}
}
\examples{
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

plot_ClaDS_phylo(tree,speciation_rates, lwd = 4, log = FALSE)

}
