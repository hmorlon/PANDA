likelihood_sgd <- function(phylo, tot_time, b, d, nu, f)
{
	# The function that everybody should use : transforms the object of type "phylo" into lambert representation, then returns the total likelihood of the tree
	lambert = Phylo2Lambert(phylo)
	lambert[1] <- tot_time
	likeli <- LikelihoodSGDFromLambert(lambert, b, d, nu, f)
	return(likeli[3]+likeli[4])
}
