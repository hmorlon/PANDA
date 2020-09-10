\name{plot_div.BipartiteEvol}
\alias{plot_div.BipartiteEvol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the output of BipartiteEvol
}
\description{
Plot the genealogies and phylogenies simulated with BipartiteEvol
}
\usage{
plot_div.BipartiteEvol(gen, spec, trait.id, lwdgen = 1, 
    lwdsp = lwdgen, scale = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{gen}{
The output of a run of make_gen.BipartiteEvol
}
  \item{spec}{
The output of a run of define_species.BipartiteEvol
}
  \item{trait.id}{
The trait dimension used to color the genealogies, phylogenies an network with trait values}
  \item{lwdgen}{
Width of the branches of the genealogies, default to 1
}
  \item{lwdsp}{
Width of the branches of the phylogenies, default to 1
}
  \item{scale}{
Optional, used to force the trait scale
}
}
\details{
The upper line shows the genealogies colored with trait values for both guilds (the number above shows the depth of the respective genealogy). 

The second line shows the phylogenies colored with trait values for both guilds (the number above shows the tip number of the respective phylogeny).

On the third line there is, from left to right, the trait distribution within individuals in guild P, trait of the individual in H as a function of the trait of the interacting individual in P, and the trait distribution within individuals in guild H (for the dimension trait.id). 

The lower line shows the quantitative interaction network, with species colored according to their mean trait value (for the dimension trait.id).
}
\references{
Maliet, O., Loeuille, N. and Morlon, H. (2020), An individual‐based model for the eco‐evolutionary emergence of bipartite interaction networks. Ecol Lett. doi:10.1111/ele.13592}
\author{
O. Maliet
}
\seealso{
\code{\link{sim.BipartiteEvol}}
}
\examples{
# run the model
set.seed(1)
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

#plot the result
plot_div.BipartiteEvol(gen,phy1, 1)
}
