\name{sim_microbiota}
\alias{sim_microbiota}

\title{
Simulate microbial alignments according to different models of Host-Microbiota evolution
}

\description{
This function simulates microbial alignments according to different models of host-microbiota evolution (strict vertical transmission, horizontal transmission, environmental acquisition).
}

\usage{
sim_microbiota(name, mu=1, n=20, seed=1,nb_cores=1,name_index,
simul,N=300,proportion_variant=0.1,path,provided_tree,...)}

\arguments{
\item{name}{
the name of the run
}

\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs) to simulate
}
\item{path}{
the path to the working directory
}

\item{nb_cores}{
a number of cores to run the analyses
}

\item{seed}{
a seed to assure the reproducibility
}

\item{simul}{
a vector of scenarii for each simulated alignment ("0" for strict vertical transmission, any positive integer for the number of host-switches transmissions, and "indep" to simulate independent evolutions and environmental acquisitions). This vector must have the same length as "name_index"
}
\item{N}{
a positive integer value representing the total number of nucleotides (variant and invariant nucleotides)
}
\item{proportion_variant}{
a value between 0 and 1 representing the proportion of variant nucleotides among N
}

\item{provided_tree}{
optional - the n-tips, binary, rooted, and ultrametric host tree to simulate the alignments on it (by default a pure-birth host tree is simulated)}

}


\value{
The function outputs microbial nucleotidic alignments according to different models of Host-Microbiota evolution. Data are well-formatted to then directly run the HOME inference (\emph{HOME_model}).

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

sim_microbiota(name="example_simulation", mu=1, n=10,
name_index=c("Simul_1","Simul_2","Simul_3"), simul=c(0,5,"indep"), 
N=300, proportion_variant=0.1)

# Run HOME (parameters estimations and model selections)

#HOME_model(name="example_simulation", name_index=c("Simul_1","Simul_2","Simul_3"), nb_tree=1000,lambda=seq(1,15), empirical=FALSE, nb_random=10)

}


