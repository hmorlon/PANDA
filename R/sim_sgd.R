sim_sgd <- function(tau, b, d, nu)
{
	lambert <- Simulate_SGD_Phylogeny(tau, 0, 2, b, d, nu)
	phylo <- read.tree(text = paste(Lambert2Newick(lambert), ";", sep=""))
	return(phylo)
}
