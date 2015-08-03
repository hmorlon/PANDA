\name{CreateGeoObject}
\alias{CreateGeoObject}

\title{
Create blank geography object
}
\description{
This function returns names of internode intervals, named descendants of each node,
and a geography object formatted in a way that can be
passed to \code{\link{fitCompModel}}

}
\usage{

CreateGeoObject(phylo,maps.object=NULL)

}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation)
}
  \item{maps.object}{
either a matrix modified from \code{phylo$edge} or a list provided from the \code{$maps} output of \code{\link{make.simmap}} (see Details and Examples)
}
}
\value{
	a list with the following components:
  \item{node.descendants}{a data frame listing node numbers, their descendants 
  (i.e., branch name), the next node, and the biogeographical region in which a node exists (if maps.object is provided)}
  \item{geography.object}{a list of matrices specifying sympatry (1) or allopatry (0) for each species pair for each internode interval (see Details)}
}

\details{
This function should be used to format the geography object so that it can be correctly
passed to the numerical integration performed in \code{\link{fitCompModel}}. 

If a \code{maps.object} is not provided, the blank \code{$geography.object} should be filled in in a block symmetrical fashion with 1s
(for sympatric species) and 0s (for allopatric species). The diagonal of each matrix
should be filled with 1s.

A \code{maps.object} can either be a matrix formed by specifying the region in which each branch specified by \code{phylo$edge} existed, or a stochastic map stored as a list in \code{$maps} from \code{\link{make.simmap}} (see Examples).
If an object from \code{\link{make.simmap}} is provided, \code{CreateGeoObject} will choose for each branch the region that was inhabited for the longest period of time along that branch.
}

\note{
The node numbers returned in \code{object$node.descendants} may not reflect those plotted using \code{nodelabels()}. The node numbers returned from
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

data(Anolis.data)
#Create a geography.object with a modified edge matrix
#First, specify which region each branch belonged to:
Anolis.regions<-c(rep("cuba",14),rep("hispaniola",17),"puerto_rico")
Anolis.map<-cbind(Anolis.data$phylo$edge,Anolis.regions)
CreateGeoObject(Anolis.data$phylo,maps.object=Anolis.map)

#Create a geography.object with a make.simmap object
#First, specify which region each branch belonged to:
geo<-c(rep("cuba",7),rep("hispaniola",9),"puerto_rico")
names(geo)<-Anolis.data$phylo$tip.label
stochastic.map<-make.simmap(Anolis.data$phylo,geo, model="ER", nsim=1)
CreateGeoObject(Anolis.data$phylo,maps.object=stochastic.map$maps)

#To create an empty geography.object:
CreateGeoObject(Anolis.data$phylo)
}
