# First attempt to improve likelihood_bd to include arguments needed for backbone anlysis

likelihood_bd_backbone <- function (phylo, tot_time, f, f.lamb.l, f.mu.l, 
          backbone, spec_times, branch_times, # arguments for backbone analysis
          cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE,
          dt = 0, cond = "crown") 
{
  if (!inherits(phylo, "phylo")) 
    stop("object \"phylo\" is not of class \"phylo\"")
#The code to compute the Likelihood
  nbtips <- Ntip(phylo)
  log_indLikelihood <- c()
  from_past <- cbind(phylo$edge, node.age(phylo)$ages)
  ages <- rbind(from_past[, 2:3], c(nbtips + 1, 0))
  ages <- ages[order(ages[, 1]), ]
  age <- max(ages[, 2])
  for (j in 1:(nbtips - 1)) {
    node <- (nbtips + j)
    edges <- phylo$edge[phylo$edge[, 1] == node, ]
    tj <- age - ages[edges[1, 1], 2]
    Psi_timevar_errap_tj <- .Psi(0, tj, f.lamb.l, f.mu.l, f, 
                                 cst.lamb = cst.lamb, cst.mu = cst.mu, expo.lamb = expo.lamb, 
                                 expo.mu = expo.mu, dt = dt)
    log_lik_tj <- log(f.lamb.l(tj)) + log(Psi_timevar_errap_tj)
    log_indLikelihood <- c(log_indLikelihood, log_lik_tj)
  }
  log_indLikelihood <- c(log_indLikelihood, log(.Psi(0, tot_time, 
                                                     f.lamb.l, f.mu.l, f, cst.lamb = cst.lamb, cst.mu = cst.mu, 
                                                     expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt)))
  
  #Type of analysis (options backbone, branch_times and spec_time should be added here)
  if (backbone == F) {
    log_data_lik <- sum(log_indLikelihood) + nbtips * log(f)
  } else if (backbone == "backbone1"){
    ########## Part similar to what might be the good way of coding for spec_times...
    #spec_lik<-1
    #for (i in 1:length(spec_times))
    #{spec_lik<-spec_lik*f.lamb.l(spec_times[i])}
    ##########
    spec_lik<-prod(f.lamb.l(spec_times))
    log_data_lik<-sum(log_indLikelihood)+nbtips*log(f)+log(spec_lik)
    
  } else if (backbone == "backbone2"){
    branch_lik<-1	
    for (k in 1:length(branch_times))
    {branch_lik<-branch_lik*f.lamb.l(branch_times[[k]][2])*.Psi(branch_times[[k]][1],branch_times[[k]][2],f.lamb=f.lamb.l,f.mu=f.mu.l,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)}
    log_data_lik<-sum(log_indLikelihood)+nbtips*log(f)+log(branch_lik)
  }
  
  #Conditioning  
  if (cond == FALSE) {
    log_final_lik <- log_data_lik
  }
  else if (cond == "stem") {
    Phi <- .Phi(tot_time, f.lamb.l, f.mu.l, f, cst.lamb = cst.lamb, 
                cst.mu = cst.mu, expo.lamb = expo.lamb, expo.mu = expo.mu, 
                dt = dt)
    log_final_lik <- log_data_lik - log(1 - Phi)
  }
  else if (cond == "crown") {
    Phi <- .Phi(tot_time, f.lamb.l, f.mu.l, f, cst.lamb = cst.lamb, 
                cst.mu = cst.mu, expo.lamb = expo.lamb, expo.mu = expo.mu, 
                dt = dt)
    log_final_lik <- log_data_lik - log(f.lamb.l(tot_time)) - 
      2 * log(1 - Phi)
  }
  return(log_final_lik)
}
