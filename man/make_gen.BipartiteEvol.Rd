\name{make_gen.BipartiteEvol}
\alias{make_gen.BipartiteEvol}
\title{
Compute the genealogies for BipartiteEvol}
\description{
Compute the genealogies from a run of BipartiteEvol
}
\usage{
make_gen.BipartiteEvol(out, treeP = NULL, treeH = NULL, verbose = TRUE)
}
\arguments{
  \item{out}{
The output of a run of sim.BipartiteEvol
}
  \item{treeP}{
Optional, a previous genealogy for clade P to which the new tree will be grafted (used if out was the continuation of a former run, see in the example)
}
  \item{treeH}{
Optional, a previous genealogy for clade H to which the new tree will be grafted (used if out was the continuation of a former run, see in the example)
}
  \item{verbose}{
Should the progression of the computation be printed?
}
}
\value{
a list object with
\item{P}{The genealogy of the clade P}
\item{H}{The genealogy of the clade H}
}
\references{
Maliet, O., Loeuille, N. and Morlon, H. (2020), An individual-based model for the eco-evolutionary emergence of bipartite interaction networks. Ecol Lett. doi:10.1111/ele.13592}
\author{
O. Maliet
}
\seealso{
\code{\link{sim.BipartiteEvol}}
}
\examples{

\dontshow{test = FALSE}
if(test){
# run the model
set.seed(1)
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