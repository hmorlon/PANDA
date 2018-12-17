\name{prepare_data_HOME}
\alias{prepare_data_HOME}

\title{
Prepare the data for running HOME 
}

\description{
This function prepares the data for the inference of HOME (format data, run the model of nucleotidic evolution...). 

}

\usage{prepare_data_HOME(iter, name, name_index, host_tree, path...)
}

\arguments{
\item{name}{
the name of the run
}
\item{iter}{
the number of iteration of a particular index from name_index
}
\item{name_index}{
a vector with the names of the different symbionts (e.g. name of the OTUs)
}
\item{host_tree}{
optional - the host tree of your study. It can also be provided in the working directory (see details).
}
\item{path}{
optional - the path to the working directory
}


}


\value{
The function creates files to run the following functions simul_bank_tree, fit_HOME...
}
\details{
Please provide the binary rooted ultrametric host tree (Newick format .tre) in the working directory. Its name must be weel-formated "host_tree_NAME.tre" (match the name of the run). A folder (“path_alignment” or "path" by default") must contain all the alignments with the filenames "alignment_NAME_OTU.fas" (FASTA format).
See reference for more details.

}
\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Modeling host-microbiota evolution for characterizing symbiont inheritance – application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{sim_microbiota}},
\code{\link{simul_bank_tree}},
\code{\link{output_results_HOME}},
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

sim_microbiota(name, name_index, simul=c(0,5,"indep"), n=10, mu=1, N=300, proportion_variant=0.1)

# Inference

# Prepare the data (format, substitution model...)
for (i in 1:3){prepare_data_HOME(iter=i,name,name_index)}

# Simulate the bank of trees
for (ksi in 1:length(seq(1,25))){simul_bank_tree(ksi,name,nb_tree=1000,lambda=seq(1,15),seed=1)}

# Infer the parameters
#for (i in 1:3){fit_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,25),nb_cores=1)}

# Plot the first outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,25),nb_tree=1000,empirical=FALSE,randomize=FALSE,raref=FALSE)}

# Perform the model selection
#for (i in 1:3){model_selection_HOME(index=name_index[i],name,nb_tree=1000,lambda=seq(1,25),nb_cores=1,seed=1)}

# Plot the final outputs
#for (i in 1:3){output_results_HOME(iter=i,name,name_index,lambda=seq(1,25),nb_tree=1000,empirical=FALSE,randomize=TRUE,raref=FALSE)}
}