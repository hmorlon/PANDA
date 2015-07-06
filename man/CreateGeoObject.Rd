\name{CreateGeoObject}
\alias{CreateGeoObject}

\title{
Create blank geography object
}
\description{
This function returns names of internode intervals, named descendants of each node,
and an unformatted geography object that are formatted in a way that can be
passed to \code{\link{fitCompModel}}

}
\usage{

CreateGeoObject(phylo)

}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
}
\value{
	a list with the following components:
  \item{node.descendants}{a data frame listing node numbers, their descendants 
  (i.e., branch name), and the next node}
  \item{extant.lineages.through.time}{a named list of extant branches that descend
  from each node (see Note)}
  \item{geography.object}{a list of empty matrices for each internode interval (see Details)}
}

\details{
This function should be used to format the geography object so that it can be correctly
passed to the numerical integration performed in \code{\link{fitCompModel}}. 

The blank \code{$geography.object} should be filled in in a block symmetrical fashion with 1s
(for sympatric species) and 0s (for allopatric species). The diagonal of each matrix
should be filled with 1s.
}

\note{
The node numbers returned in \code{object$node.descendants} and \code{object$extant.lineages.through.time}
may not reflect those plotted using \code{nodelabels()}. The node numbers returned from
this function start at \code{length(phylo$tip.label)+1} and increase chronologically.
}

\references{
	Drury, J., Clavel, J., and Morlon, H. (submitted) Estimating the effect of competition on trait evolution using maximum likelihood inference.
}

\author{Jonathan Drury \href{mailto:drury@biologie.ens.fr}{drury@biologie.ens.fr}}

\seealso{
\code{\link{fitCompModel}}
}

\examples{

#To create a geography.object
tree<-rcoal(5)
tree_internal.nodes<-CreateGeoObject(tree)
geography.matrix<-tree_internal.nodes$geography.object
geography.matrix[[1]][,1]<-c(1,0)
geography.matrix[[1]][,2]<-c(0,1)
geography.matrix[[2]][,1]<-c(1,1,0)
geography.matrix[[2]][,2]<-c(1,1,0)
geography.matrix[[2]][,3]<-c(0,0,1)
geography.matrix[[3]][,1]<-c(1,1,0,0)
geography.matrix[[3]][,2]<-c(1,1,0,0)
geography.matrix[[3]][,3]<-c(0,0,1,1)
geography.matrix[[3]][,4]<-c(0,0,1,1)
geography.matrix[[4]][,1]<-c(1,1,1,0,0)
geography.matrix[[4]][,2]<-c(1,1,1,0,0)
geography.matrix[[4]][,3]<-c(1,1,1,0,0)
geography.matrix[[4]][,4]<-c(0,0,0,1,1)
geography.matrix[[4]][,5]<-c(0,0,0,1,1)
}
