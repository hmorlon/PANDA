fit_sgd <- function(phylo, tot_time, par_init, f=1)
{
	#par_init contains the initial values of parameters (b, d, nu), phylo is an object of class phylo
	lambert = Phylo2Lambert(phylo)
	lambert[1] <- tot_time

	optim_LikelihoodSGD <- function(par)
	{
		likeli = LikelihoodSGDFromLambert(lambert, abs(par[1]), abs(par[2]), abs(par[3]), f)
		L = likeli[3]+likeli[4]
		return(-log(L))
	}
	
	inference <- optim(par_init, optim_LikelihoodSGD)$par
	return( c(birth = abs(inference[1]), growth = abs(inference[1])-abs(inference[2]), mutation = abs(inference[3])) )
}
