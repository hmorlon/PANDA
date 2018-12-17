\name{HOME_model}
\alias{HOME_model}

\title{
HOst-Microbiota Evolution model
}

\description{
This function runs HOME for all the different OTUs (called "index"") on the host tree.

}

\usage{
HOME_model(name, name_index, provided_tree, nb_tree=5000, lambda=seq(1,25),
raref=FALSE, empirical=TRUE, randomize=TRUE, nb_random=10, path, path_alignment, seed=3, nb_cores=1)
}

\arguments{
\item{name}{
the name of the run
}
\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs)
}

\item{provided_tree}{
optional - the host tree of your study. It can also be provided as saved in the working directory (see details).
}

\item{nb_tree}{
a number of tree for Monte-Carlo estimation of the number of switches (a low value will give inaccurate results whereas a high value will increase the computation time)
}
\item{lambda}{
a vector of integer values of number of switches to test during estimations
}
\item{raref}{
if TRUE rarefactions on the number of trees are performed (i.e. to test if nb_tree is large enough)
}
\item{empirical}{
a boolean indicating if the run is done on empirical data or on simulations (in the latter case, the output will be different)
}
\item{randomize}{
a boolean indicating to perform or not the model selection testing independent evolutions (it takes time)
}

\item{nb_random}{
a integer value of the number of randomizations (R) to perform during the model selection testing independent evolutions
}

\item{path}{
optional - the path to the working directory
}
\item{path_alignment}{
optional - the path toward the folder containing the different OTU alignments (by default, the working directory)
}
\item{seed}{
a seed to assure the reproducibility
}
\item{nb_cores}{
a number of cores to run the analyses (ideally, it should be equal to the length of lambda for an optimal speed)
}
}


\value{
The function gives all the results of HOME in separate files and figures (see inside folders results/ and figures/) and sum up all the results for each OTU in a HTML file (see in folder figures/). The function HOME_model contains the functions \emph{prepare_data_HOME}, \emph{simul_bank_tree}, \emph{fit_HOME}, \emph{model_selection_HOME}, and \emph{output_results_HOME} that can be run succesively.
}
\details{

Alignments must be provided in FASTA format and their name must be formatted "alignment_NAME_INDEX.fas" (match the name of the run and the name of the index of the alignment).

If you do not provide a host tree directly in the function, please provide the binary rooted ultrametric host tree (Newick format .tre) in the working directory. The name of the file name must be well-formated "host_tree_NAME.tre" (match the name of the run).

See reference for more details.

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Characterizing symbiont inheritance during host-microbiota evolution: application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{sim_microbiota}}
}
\examples{
# Some examples may take a little bit of time. Be patient!

### ON SIMULATIONS 

# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)

sim_microbiota(name="example_simulation", name_index=c("Simul_1","Simul_2","Simul_3"),
simul=c(0,5,"indep"), n=10, mu=1, 
N=300, proportion_variant=0.1)

# Run HOME (parameters estimations and model selections)

#HOME_model(name="example_simulation", name_index=c("Simul_1","Simul_2","Simul_3"), nb_tree=1000, lambda=seq(1,15), empirical=FALSE, nb_random=10)


### EMPIRICAL APPLICATION: example on the great apes microbiota

# Load the data
example_great_apes_microbiota(name="great_apes")

# Run HOME
#HOME_model(name="great_apes", name_index=c("OTU0001","OTU0002","OTU0003"), nb_tree=500, lambda=seq(1,25), nb_random=10)

# NB: this should be run on a larger number of tree!


}





























