\name{sim.BipartiteEvol}
\alias{sim.BipartiteEvol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the BipartiteEvol model}
\description{
Simulateof the BipartiteEvol model from Maliet et al. (2020)
}
\usage{
sim.BipartiteEvol(nx, ny = nx, NG, dSpace = Inf, D = 1, muP,
muH, alphaP = 0, alphaH = 0, iniP = 0, iniH = 0, nP = 1, nH = 1, 
rP = 1, rH = 1, effect = 1, verbose = 100, thin = 1, P = NULL, H = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nx}{
Size of the grid (the grid has size nx * ny)
}
  \item{ny}{
Size of the grid (default to nx, the grid has size nx * ny)
}
  \item{NG}{
Number of time step the model is run
}
  \item{dSpace}{
Size of the dispersal kernel (default to Inf, meaning there are no restrictions on dispersion)}
  \item{D}{
Dimention of the trait space (default to 3)}
  \item{muP}{
Mutation probability at reproduction for the individuals of clade P}
  \item{muH}{
Mutation probability at reproduction for the individuals of clade H}
  \item{alphaP}{
alpha parameter for clade P (1/alpha is the niche width)}
  \item{alphaH}{
alpha parameter for clade H (1/alpha is the niche width)}
  \item{iniP}{
Initial trait value for the individuals in clade P}
  \item{iniH}{
Initial trait value for the individuals in clade P}
  \item{nP}{
Number of individuals of clade P killed at each time step
}
  \item{nH}{
Number of individuals of clade H killed at each time step
}
  \item{rP}{
r parameter for clade P (r is the ratio between the fitness maximum and minimum)}
  \item{rH}{
r parameter for clade H (r is the ratio between the fitness maximum and minimum)
}
  \item{effect}{
Standard deviation of the trait mutation kernel
}
  \item{verbose}{
The simulation }
  \item{thin}{
The number of iterations between two recording of the state of the model (default to 1)
}
  \item{P}{
Optionnal, used to continue one precedent run: traits of the individuals of clade P at the end of the precedent run
}
  \item{H}{
Optionnal, used to continue one precedent run: traits of the individuals of clade H at the end of the precedent run
}
}

\value{
a list with 
\item{Pgenealogy}{The genalogy of clade P}
\item{Hgenealogy}{The genalogy of clade H}
\item{xP}{The trait values at each time step for clade P}
\item{xH}{The trait values at each time step for cladeH}
\item{P}{The trait values at present for clade P}
\item{H}{The trait values at present for clade P}
\item{Pmut}{The number of new mutations at each time step for clade P}
\item{Hmut}{The number of new mutations at each time step for clade H}
\item{iniP}{The initial trait values for the individuals of clade P used in the simulation}
\item{iniH}{The initial trait values for the individuals of clade H used in the simulation}
\item{thin.factor}{The thin value used in the simulation}
}
\references{
Maliet, O., Loeuille, N. and Morlon, H. (2020), An individual‐based model for the eco‐evolutionary emergence of bipartite interaction networks. Ecol Lett. doi:10.1111/ele.13592}
\author{
O. Maliet
}
\examples{
# run the model
set.seed(1)
mod = sim.BipartiteEvol(nx = 8,ny = 4,NG = 500,
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
plot_net.BipartiteEvol(gen,phy1,trait.id, net,mod, nx = 10, spatial = FALSE)


## add time steps to a former run
seed=as.integer(10)
set.seed(seed)

mod = sim.BipartiteEvol(nx = 8,ny = 4,NG = 500,
                        D = 3, muP = 0.1 , muH = 0.1,
                        alphaP = 0.12,alphaH = 0.12,
                        rP = 10, rH = 10,
                        verbose = 100, thin = 5,
                        P=mod$P,H=mod$H)  # former ru output

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

