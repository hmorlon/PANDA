\name{simul_bank_tree}
\alias{simul_bank_tree}

\title{
Simulate a bank of trees for running HOME.
}

\description{
This function simulates the bank of symbiont trees, used for running HOME, for a given number of host-switches.

}

\usage{simul_bank_tree(nb_ksi,name,
provided_tree=NULL,nb_tree=10000,
lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),
seed=1)}

\arguments{
\item{nb_ksi}{
a positive integer which represents the given number of host-switches to simulate on the host tree
}

\item{name}{
the name of the run
}

\item{provided_tree}{
optional - the host tree of your study. It can also be provided as saved in the working directory}

\item{nb_tree}{
a number of trees to simulate
}
\item{lambda}{
a vector of integer values of number of switches to test during estimations
}
\item{seed}{
a seed to assure the reproducibility
}

}


\value{
The function outputs trees in the folder simulated_trees/ for each value of number of switches. Then, use the function \emph{fit_HOME} to infer the results.
}
\details{

The function \emph{prepare_data_HOME} must be run before.
See reference for more details.

}
\references{
Perez-Lamarque B, Morlon H (2019). Characterizing symbiont inheritance during host-microbiota evolution: Application to the great apes gut microbiota. Molecular Ecology Resources 19:1659-1671.
}
\author{Benoit Perez-Lamarque}
\seealso{
\code{\link{sim_microbiota}},
\code{\link{prepare_data_HOME}},
\code{\link{output_results_HOME}},
\code{\link{model_selection_HOME}},
\code{\link{fit_HOME}},
\code{\link{HOME_model}}
}
\examples{
# Some examples may take a little bit of time. Be patient!

# Simulate 3 microbial alignments on a host tree
# (1 is vertically transmitted, 1 is transmitted with 5 host-switches,
#  and 1 is environmentally acquired)

name="example_simulation"
name_index=c("Simul_1","Simul_2","Simul_3")
path=getwd()

#sim_microbiota(name, name_index, simul=c(0,5,"indep"), n=10, mu=1, N=300, proportion_variant=0.1)

# Inference

# Prepare the data (format, substitution model...)
#for (i in 1:3){prepare_data_HOME(iter=i,name,name_index)}

# Simulate the bank of trees
#for (ksi in 1:length(seq(1,25))){simul_bank_tree(ksi,name,nb_tree=1000,
#lambda=seq(1,15),seed=1)}

# Infer the parameters
#for (i in 1:3){fit_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,25),nb_cores=1)}

# Plot the first outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,25),
#nb_tree=1000,empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
#for (i in 1:3){model_selection_HOME(index=name_index[i],name,nb_tree=1000,
#lambda=seq(1,25),nb_cores=1,seed=1)}

# Plot the final outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,25),nb_tree=1000,
#empirical=FALSE,randomize=TRUE,raref=FALSE)}
}
