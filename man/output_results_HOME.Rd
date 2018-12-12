\name{output_results_HOME}
\alias{output_results_HOME}

\title{
Representing all the results of HOME model
}

\description{
This function outputs all the results into a HTML file with plots and tables of the HOME model. The results are reported in separate file for each OTU (index).

}

\usage{output_results_HOME(iter,name,name_index,lambda,nb_tree=1000,empirical,randomize,raref,...)}

\arguments{
\item{name}{
the name of the run
}
\item{iter}{
the number of iteration of a particular index from name_index
}
\item{nb_cores}{
a number of cores to run the analyses
}
\item{nb_tree}{
a number of tree for Monte-Carlo estimation of the number of switches (a low value will give inaccurate results whereas a high value will increase the computation time)
}
\item{lambda}{
a vector integer values of number of switches to test during estimations
}
\item{raref}{
if TRUE rarefactions on the number of trees have been performed (i.e. to test if nb_tree is large enough)
}
\item{empirical}{
a boolean indicating if the run has been done on empirical data or on simulations (in the latter case, the output will be different)
}
\item{randomize}{
a boolean indicating wheter the model selection testing independent evolutions has been performed yet.
}

}


\value{
This function sums up all the results for each index in a HTML file (see in folder figures/).}
\details{

The functions \emph{prepare_data_HOME}, \emph{simul_bank_tree}, and \emph{fit_HOME} must be run before.
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
\code{\link{model_selection_HOME}},
\code{\link{fit_HOME}},
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





























