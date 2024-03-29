\name{define_species.BipartiteEvol}
\alias{define_species.BipartiteEvol}
\title{
Build the phylogenies for BipartiteEvol}
\description{
Build the phylogenies from the output of BipartiteEvol and the corresponding genealogies}
\usage{
define_species.BipartiteEvol(genealogy, threshold = 1, 
      distanceH = NULL, distanceP = NULL, verbose = T,
      monophyly = TRUE, seed = NULL)
}
\arguments{
  \item{genealogy}{
The output of a run of make_gen.BipartiteEvol
}
  \item{threshold}{
The species definition ratchet (s)
}
  \item{distanceH}{
Distance (ie nb of mutations) matrix between the individual of clade H
}
  \item{distanceP}{
Distance (ie nb of mutations) matrix between the individual of clade P}
  \item{verbose}{
Should the progression of the computation be printed?
}
  \item{monophyly}{ 
Should the species delineations be strictly monophyletic species (TRUE - default) or not (FALSE)? If not, the threshold must be equal to 1.
}
  \item{seed}{
If monophyly==FALSE, the seed is used to pick one representative individual per (potentially non-monophyletic) species. 
}
}
\value{
a list with
\item{P}{The species identity of each individual in guild P}
\item{H}{The species identity of each individual in guild H}
\item{Pphylo}{The phylogeny for guild P}
\item{Hphylo}{The phylogeny for guild H}
}

\details{
If monophyly==TRUE, species delineation is performed using the model of Speciation by Genetic Differentiation (Manceau et al., 2015) where the `threshold` (the number of mutations needed to belong to different species) can vary. It results in monophyletic species.
If monophyly==FALSE, we consider that each new mutation (i.e. each new combination of traits) gives rise to a new species (Perez-Lamarque et al., 2021). As a result, species are not necessarily formed by a monophyletic group of individuals. 
}

\references{

Manceau, M., A. Lambert, and H. Morlon. (2015). Phylogenies support out-of-equilibrium models of biodiversity. Ecology letters 18:347–356.

Maliet, O., Loeuille, N. and Morlon, H. (2020). An individual-based model for the eco-evolutionary emergence of bipartite interaction networks. Ecol Lett. doi:10.1111/ele.13592

Perez‐Lamarque, B., Maliet, O., Pichon B., Selosse, M-A., Martos, F., Morlon H. (2021). Do closely related species interact with similar partners? Testing for phylogenetic signal in bipartite interaction networks. bioRxiv. doi: https://doi.org/10.1101/2021.08.30.458192

}
\author{
O. Maliet & B. Perez-Lamarque
}
\seealso{
\code{\link{sim.BipartiteEvol}}
}
\examples{
# run the model
set.seed(1)

\dontshow{test <- FALSE}
if(test){

mod = sim.BipartiteEvol(nx = 8,ny = 4,NG = 800,
                        D = 3, muP = 0.1 , muH = 0.1,
                        alphaP = 0.12,alphaH = 0.12,
                        rP = 10, rH = 10,
                        verbose = 100, thin = 5)

#build the genealogies
gen = make_gen.BipartiteEvol(mod)
plot(gen$H)

#compute the phylogenies
phy1 = define_species.BipartiteEvol(gen,threshold=1)

#plot the result
plot_div.BipartiteEvol(gen,phy1, 1)

#build the network
net = build_network.BipartiteEvol(gen, phy1)

trait.id = 1
plot_net.BipartiteEvol(gen,phy1,trait.id, net,mod, nx = nx, spatial = FALSE)


## add time steps to a former run
seed=as.integer(10)
set.seed(seed)

mod = sim.BipartiteEvol(nx = 8,ny = 4,NG = 200,
                        D = 3, muP = 0.1 , muH = 0.1,
                        alphaP = 0.12,alphaH = 0.12,
                        rP = 10, rH = 10,
                        verbose = 100, thin = 5,
                        P=mod$P,H=mod$H)  # former run output

# update the genealogy
gen = make_gen.BipartiteEvol(mod,
                             treeP=gen$P, treeH=gen$H)

# update the phylogenies...
phy1 = define_species.BipartiteEvol(gen,threshold=1)

#... and the network
net = build_network.BipartiteEvol(gen, phy1)

trait.id = 1
plot_net.BipartiteEvol(gen,phy1,trait.id, net,mod, nx = 10, spatial = FALSE)

}

}

