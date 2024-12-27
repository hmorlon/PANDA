\name{plot_net.BipartiteEvol}
\alias{plot_net.BipartiteEvol}
\title{
Plot the output of BipartiteEvol
}
\description{
Plot the genealogies, phylogenies and interaction network simulated with BipartiteEvol
}
\usage{
plot_net.BipartiteEvol(gen, spec, trait.id, link, 
    out, lwdgen = 1, lwdsp = lwdgen, scale = NULL, 
    nx = NULL, cor = FALSE, network.method = "bipartite", 
    spatial = FALSE)
}
\arguments{
  \item{gen}{
The output of a run of make_gen.BipartiteEvol
}
  \item{spec}{
The output of a run of define_species.BipartiteEvol
}
  \item{trait.id}{
The trait dimension used to color the genealogies, phylogenies an network with trait values}
  \item{out}{
The output of a run of sim.BipartiteEvol
}
\item{link}{
The output of a run of sim.BipartiteEvol
}
  \item{lwdgen}{
Width of the branches of the genealogies, default to 1
}
  \item{lwdsp}{
Width of the branches of the phylogenies, default to 1
}
  \item{scale}{
Optional, used to force the trait scale
}
  \item{nx}{
Grid size parameter used in sim.BipartiteEvol. If NULL, squrt(N) is used, where N is the number of individuals in a guild}
  \item{cor}{
If FALSE (the default), the middle panel displays the interraction network with species positionned in trait space. If TRUE, it shows all the individual in trait space}
  \item{network.method}{
How should the network be plotted? Can be "bipartite" (the default) or "matrix"}
  \item{spatial}{
Should the grid with trait values of the individual of both guilds been shown? Default to FALSE}
}

\value{Plot outputs of BipartiteEvol} 

\details{
The upper line shows the genealogies colored with trait values for both guilds (the number above shows the depth of the respective genealogy). 

The second line shows the phylogenies colored with trait values for both guilds (the number above shows the tip number of the respective phylogeny).

On the third line there is, from left to right, the trait distribution within individuals in guild P (for the dimension trait.id), the interraction network with species positionned in trait space (if cor = T), and the trait distribution within individuals in guild H (for the dimension trait.id). 

The lower line shows the quantitative interaction network, with species colored according to their mean trait value (for the dimension trait.id).
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
# run the model
set.seed(1)

\dontshow{test = FALSE}
if(test){
mod = sim.BipartiteEvol(nx = 8,ny = 4,NG = 1000,
                        D = 3, muP = 0.1 , muH = 0.1,
                        alphaP = 0.12,alphaH = 0.12,
                        rP = 10, rH = 10,
                        verbose = 100, thin = 5)

#build the genealogies
gen = make_gen.BipartiteEvol(mod)
plot(gen$H)

#compute the phylogenies
phy1 = define_species.BipartiteEvol(gen,threshold=1)

#build the network
net = build_network.BipartiteEvol(gen, phy1)

trait.id = 1
plot_net.BipartiteEvol(gen,phy1,trait.id, net,mod, nx = 10, spatial = FALSE)
}
}
