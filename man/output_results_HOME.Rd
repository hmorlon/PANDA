\name{output_results_HOME}
\alias{output_results_HOME}

\title{
Representing all the results of HOME
}

\description{
This function outputs all the results of HOME into a HTML file with figures and tables. The results are reported in one separate file for each OTU ("index").

}

\usage{output_results_HOME(iter,name,name_index,lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),
nb_tree=10000,empirical,randomize,raref,nb_random=10,figure=FALSE,...)}

\arguments{
\item{iter}{
the number of iteration of a particular index from name_index
}
\item{name}{
the name of the run
}
\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs)
}
\item{nb_tree}{
a number of tree for Monte Carlo estimation of the number of switches (a low value will give inaccurate results whereas a high value will increase the computation time)
}
\item{lambda}{
a vector integer values of number of switches to test during estimations
}

\item{empirical}{
a boolean indicating if the run has been done on empirical data or on simulations (in the latter case, the output will be different)
}
\item{randomize}{
a boolean indicating wheter the model selection testing independent evolutions has been performed yet.
}
\item{raref}{
if TRUE rarefactions on the number of trees have been performed (i.e. to test if nb_tree is large enough) and a separate figure with rarefactions curves will be stored in the folder "figures/"
}

\item{nb_random}{
a integer value of the number of randomizations (R) performed during the model selection testing independent evolutions.
}

\item{figure}{
indicates whether it automatically outputs figures (default: FALSE).
}

\item{...}{
optional - other arguments to be passed.
}

}


\value{
This function sums up all the results for each index in a HTML file (see in folder "figures/").}
\details{

The functions \emph{prepare_data_HOME}, \emph{simul_bank_tree}, and \emph{fit_HOME} must be run before.
See reference for more details.

}
\references{
Perez-Lamarque B, Morlon H (2019). Characterizing symbiont inheritance during host-microbiota evolution: Application to the great apes gut microbiota. Molecular Ecology Resources 19:1659-1671.
}
\author{Benoit Perez-Lamarque}
\seealso{
\code{\link{sim_microbiota}},
\code{\link{prepare_data_HOME}},
\code{\link{simul_bank_tree}},
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
#lambda=seq(1,25),seed=1)}

# Infer the parameters
#for (i in 1:3){fit_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,25),nb_cores=1)}

# Plot the first outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,25),nb_tree=1000,
#empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
#for (i in 1:3){model_selection_HOME(index=name_index[i],name,nb_tree=1000,
#lambda=seq(1,25),nb_cores=1,seed=1)}

# Plot the final outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,25),nb_tree=1000,
#empirical=FALSE,randomize=TRUE,raref=FALSE)}
}
