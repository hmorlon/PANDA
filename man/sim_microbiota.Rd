\name{sim_microbiota}
\alias{sim_microbiota}

\title{
Simulate microbial alignments according to different models of Host-Microbiota evolution
}

\description{
This function simulates microbial alignments according to different models of host-microbiota evolution (strict vertical transmission, horizontal transmission, environmental acquisition).
}

\usage{
sim_microbiota(name, name_index, simul, provided_tree, n=20, mu=1, N=300, proportion_variant=0.1, seed=1, nb_cores=1, path,...)

sim_microbiota(name, name_index, simul, provided_tree=NULL, n=20, mu=1, N=300, proportion_variant=0.1, model="uniform", mean=0.5, sd=0.01, host_signal=10, geo_signal=0, stochastic_map=NULL, path=getwd(), seed=1, nb_cores=1,...)

}

\arguments{
\item{name}{
the name of the run
}

\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs) to simulate
}

\item{simul}{
a vector of scenarii for each simulated alignment ("0" for strict vertical transmission, any positive integer for the number of host-switches transmissions, and "indep" to simulate independent evolutions and environmental acquisitions). This vector must have the same length as "name_index"
}

\item{provided_tree}{
optional - the n-tips, binary, rooted, and ultrametric host tree to simulate the alignments on it (by default a pure-birth host tree is simulated)}

\item{n}{
optional - a positive integer value representing the number of tips of the host tree to simulate (need to be provided if you do not already provide a host tree)
}

\item{mu}{
a positive value representing the simulated substitution rate. 
}

\item{N}{
a positive integer value representing the total number of nucleotides (variable and unvariable nucleotides)
}

\item{proportion_variant}{
a value between 0 and 1 representing the proportion of variable nucleotides among N
}

\item{model}{
Several ways to simulate host-switches on the host phylogeny are  available: "model=uniform" simulates host-switches uniformly on the host phylogeny, "model=temporal" simulates host-switches preferentially at a given time in the past, "model=host_dependent" simulates host-switches preferentially between closely related hosts lineages, and "model=geo_dependent" simulates host-switches preferentially between host lineages sharing the same geographic area. (default="uniform").
}

\item{mean}{
optional – for "model=temporal", the host-switches are simulated at a particular time in the past. The temporal distribution of host-switches follows a normal distribution centered on the "mean" and with a standard deviation "sd". "mean" must be a number between 0 and 1 indicating at what time in the past host-switches are more likely: "means=0" simulates host-switches closed to the root of the host phylogeny, whereas "means=1" simulates host-switches closed to the tips.
}

\item{sd}{
optional – for "model=temporal", the host-switches are simulated at a particular time in the past. The temporal distribution of host-switches follows a normal distribution centered on the "mean" and with a standard deviation "sd" ("sd" must be positive).
}

\item{host_signal}{
optional – for "model= host_dependent", the host-switches are preferentially simulated between host lineages that are closely related. "host_signal" must be a positive number. If "host_signal=0", switches are uniformly distributed, whereas the phylogenetically preference increases with the "host_signal" value.
}

\item{geo_signal}{
optional – for "model=geo_dependent", the host-switches are preferentially simulated between host lineages sharing the same geographic area. "geo_signal" must be between 0 and 1, it indicates the probability for a switch to occur between areas divided by the probability for a switch to occur within an area. If "geo_signal"=1, switches are uniformly distributed, whereas if "geo_signal"=0, switches can only happen between host lineages sharing the same geographic area. 
}

\item{stochastic_map}{
optional – for "model=geo_dependent", the host-switches are preferentially simulated between host lineages sharing the same geographic area. "stochastic_map" represents a stochastic mapping that indicates the ancestral biogeography of the host lineages (it indicates which host lineages occupied which geographic areas at a given time, and thus indicates which lineages are sharing the same vs. different areas). "stochastic_map" must be an A object of class "simmap" (R-package phytools, Revell, 2012).
}

\item{seed}{
a seed to assure the reproducibility
}

\item{nb_cores}{
a number of cores to run the analyses
}


\item{path}{
optional - the path to the working directory
}

}


\value{
The function outputs microbial nucleotidic alignments according to different models of Host-Microbiota evolution. Data are well-formatted to then directly run HOME (\emph{HOME_model}).

}
\details{

If you provide a binary rooted ultrametric host tree, it must be in Newick format (.tre).

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Characterizing symbiont inheritance during host-microbiota evolution : application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{HOME_model}}
}
\examples{

# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)

sim_microbiota(name="example_simulation", name_index=c("Simul_1","Simul_2","Simul_3"),
simul=c(0,5,"indep"), n=10, mu=1, 
N=300, proportion_variant=0.1)

# Run HOME (parameters estimations and model selections)

#HOME_model(name="example_simulation", name_index=c("Simul_1","Simul_2","Simul_3"), nb_tree=1000, lambda=seq(1,25), empirical=FALSE, nb_random=10)


}


