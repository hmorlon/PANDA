\name{CreateGeoObject_BioGeoBEARS}
\alias{CreateGeoObject_BioGeoBEARS}

\title{
Create biogeography object using a stochastic map from BioGeoBEARS
}
\description{
Create biogeography object using a stochastic map created from any model in BioGeoBEARS (see documentation in BioGeoBEARS package).

}
\usage{
CreateGeoObject_BioGeoBEARS( full.phylo, trimmed.phylo = NULL, ana.events,
clado.events, stratified=FALSE, simmap.out=FALSE)
}


\arguments{
  \item{full.phylo}{
the object of type 'phylo' (see ape documentation) that was used to construct the stochastic map in BioGeoBEARS
}
  \item{trimmed.phylo}{
if the desired biogeography object excludes some species that were initially included in the stochastic map, this specifies a phylo object for the trimmed set of species
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
  \item{stratified}{
logical indicating whether output should be a stochastic map (simmap) object (see note)
}

}
\value{
	a list with the following components:
  \item{geography.object}{a list of matrices specifying sympatry (1) or allopatry (0) for each species pair for each internode interval (see Details)}
  \item{times}{a vector containing the time since the root of the tree at which nodes or changes in biogeography occur (used internally in other functions)}
  \item{spans}{a vector specifying the distances between times (used internally in other functions)}
}

\details{

Note: generating a stochastic map output using \code{simmap.out=TRUE} and passing to \link{fit_t_comp} for diversity dependent models with biogeography greatly speeds up model 
fitting compared to output generated when \code{simmap.out=FALSE}. This cannot be used for matching competition or any two-regime models with biogeography.

}


\references{
	Drury, J., Clavel, J., Manceau, M., and Morlon, H. 2016. Estimating the effect of competition on trait evolution using maximum likelihood inference. \emph{Systematic Biology} doi 10.1093/sysbio/syw020

	Matzke, N. 2014. Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. \emph{Systematic Biology} 63:951-970.
}

\author{

	Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

}

\seealso{
\code{\link{fit_t_comp}}
\code{\link{CreateGeoObject}}
}

\examples{

data(BGB.examples)

\donttest{


##Example with a non-stratified tree

Canidae.geography.object<-CreateGeoObject_BioGeoBEARS(full.phylo=BGB.examples$Canidae.phylo,
ana.events=BGB.examples$Canidae.ana.events, clado.events=BGB.examples$Canidae.clado.events)

#on a subclade
Canidae.trimmed<-drop.tip(BGB.examples$Canidae.phylo 
							,BGB.examples$Canidae.phylo$tip.label[1:9])
							
Canidae.trimmed.geography.object<-CreateGeoObject_BioGeoBEARS(
full.phylo=BGB.examples$Canidae.phylo, trimmed.phylo=Canidae.trimmed, 
ana.events=BGB.examples$Canidae.ana.events, clado.events=BGB.examples$Canidae.clado.events)

##Example with a stratified tree

Ochotonidae.geography.object<-CreateGeoObject_BioGeoBEARS( 
full.phylo = BGB.examples$Ochotonidae.phylo, ana.events = BGB.examples$Ochotonidae.ana.events,
clado.events = BGB.examples$Ochotonidae.clado.events, stratified = TRUE)

#on a subclade
Ochotonidae.trimmed<-drop.tip(BGB.examples$Ochotonidae.phylo, 
BGB.examples$Ochotonidae.phylo$tip.label[1:9])
								
Ochotonidae.trimmed.geography.object<-CreateGeoObject_BioGeoBEARS(
full.phylo=BGB.examples$Ochotonidae.phylo, trimmed.phylo=Ochotonidae.trimmed, 
ana.events=BGB.examples$Ochotonidae.ana.events, 
clado.events=BGB.examples$Ochotonidae.clado.events, stratified=TRUE)

}

}
