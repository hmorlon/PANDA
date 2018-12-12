\name{fit_HOME}
\alias{fit_HOME}

\title{
General inference of HOME 
}

\description{
This function infers parameters (substitution rate and number of host-switches) of the HOME model for one given OTU ("index").

}

\usage{fit_HOME(index,name,nb_tree=1000,lambda,
nb_cores=1,raref,...)}

\arguments{
\item{name}{
the name of the run
}
\item{index}{
the name of the particular index from name_index (e.g. name of one OTU)
}
\item{nb_cores}{
a number of cores to run the analyses (ideally, it should be equal to the length of lambda for an optimal speed)
}
\item{nb_tree}{
a number of tree for Monte-Carlo estimation of the number of switches (a low value will give inaccurate results whereas a high value will increase the computation time)
}
\item{lambda}{
a vector of integer values of number of switches to test during estimations
}
}


\value{
The function infers the parameters (subtitution rate and number of switches) for one given index. Use the function output_results to interpret them.
}
\details{

The functions prepare_data and simul_bank_tree must be run before.
See reference for more details.

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Characterizing symbiont inheritance during host-microbiota evolution : application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{sim_microbiota}},
\code{\link{prepare_data_HOME}},
\code{\link{simul_bank_tree}},
\code{\link{output_results_HOME}},
\code{\link{model_selection_HOME}},
\code{\link{HOME_model}}
}
\examples{
# Some examples may take a little bit of time. Be patient!


# Simulate 3 microbial alignments on a host tree (1 is vertically transmitted, 1 is transmitted with 5 host-switches and 1 is environmentally acquired)

name="example_simulation"
name_index=c("Simul_1","Simul_2","Simul_3")
path=getwd()

sim_microbiota(name, mu=1, n=10,name_index, simul=c(0,5,"indep"), N=300, proportion_variant=0.1)

# Inference

# Prepare the data (format, substitution model...)
for (i in 1:3){prepare_data_HOME(iter=i,name,name_index,path)}

# Simulate the bank of trees
for (ksi in 1:length(seq(1,15))){simul_bank_tree(ksi,name,nb_tree=1000,lambda=seq(1,15),seed=1)}

# Infer the parameters
#for (i in 1:3){fit_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1)}

# Plot the first outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
#for (i in 1:3){model_selection_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,15),nb_cores=1,seed=1)}

# Plot the final outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,15),nb_tree=1000,empirical=FALSE,randomize=TRUE,raref=FALSE)}

}





























