\name{model_selection_HOME}
\alias{model_selection_HOME}

\title{
Perform the model selections of HOME
}

\description{
This function performs the model selections of HOME testing the null hypotheses of strict vertical transmission and independent evolution.

}

\usage{model_selection_HOME(index,name,nb_tree=10000,
lambda=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25),
nb_cores=1,seed=1,nb_random=10,tolerance=0.05,overwrite=TRUE,...)
}

\arguments{
\item{index}{
the name of the particular index from name_index (e.g. name of one OTU)
}

\item{name}{
the name of the run
}

\item{nb_tree}{
a number of tree for Monte Carlo estimation of the number of switches (a low value will give inaccurate results whereas a high value will increase the computation time)
}

\item{nb_random}{
a integer value of the number of randomizations (R) to perform during the model selection testing independent evolutions
}

\item{lambda}{
a vector of integer values of number of switches to test during estimations
}
\item{seed}{
a seed to assure the reproducibility
}

\item{nb_cores}{
a number of cores to run the analyses (ideally, it should be equal to the length of lambda for an optimal speed)
}
\item{tol}{
the desired accuracy of the \emph{optimize} function (to estimate mu). A low value will give a more accurate estimate, but will take more computational time.
}
\item{overwrite}{
indicates whether the any previous results should be overwritten by the current run (default: TRUE).
}
\item{...}{
optional - other arguments to be passed.
}

}


\value{
The function outputs different files. Use the function \emph{output_results_HOME} to interpret them.}
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
\code{\link{output_results_HOME}},
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
#for (ksi in 1:length(seq(1,25))){simul_bank_tree(ksi,name,nb_tree=1000,lambda=seq(1,25),seed=1)}

# Infer the parameters
#for (i in 1:3){fit_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,25),
#nb_cores=1)}

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

