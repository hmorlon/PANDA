likelihood_bd <- function(phylo,tot_time,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=1e-3,cond="crown")
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  nbtips <- Ntip(phylo)
  log_indLikelihood <- c()
  from_past <- cbind(phylo$edge,node.age(phylo)$ages)
  ages <- rbind(from_past[,2:3],c(nbtips+1,0))
  ages <- ages[order(ages[,1]),]
  age <- max(ages[,2])
  # the objective is to define once the r.int & r.int.int functions
  psi_phi_elements <- .defineRintRintint(f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,dt, tot_time)
  for (j in 1:(nbtips-1))
  {
    node <- (nbtips+j)
    edges <- phylo$edge[phylo$edge[,1]==node,]
    tj <- age-ages[edges[1,1],2]
    rst <- psi_phi_elements$r.int(0,tj)
    rist <- psi_phi_elements$r.int.int(0,tj)
    ri0s <- psi_phi_elements$r.int.int(0,0)
    Psi <- rst - 2*log(abs(1+rist/(1/f+ri0s)))
    log_lik_tj <- log(f.lamb(tj)) + Psi
    log_indLikelihood <- c(log_indLikelihood, log_lik_tj)
  }
  rst <- psi_phi_elements$r.int(0,tot_time)
  rist <- psi_phi_elements$r.int.int(0,tot_time)
  ri0s <- psi_phi_elements$r.int.int(0,0)
  Psi <- rst - 2*log(abs(1+rist/(1/f+ri0s)))
  log_lik_tot_time <- log(f.lamb(tj)) + Psi
  log_indLikelihood <- c(log_indLikelihood, log_lik_tot_time)
  # Compute total likelihood
  log_data_lik <- sum(log_indLikelihood)+nbtips*log(f)
  if (cond==FALSE)
  {
    log_final_lik <- log_data_lik
  }

  else if (cond=="stem")
  {
    rit <- psi_phi_elements$r.int(0,tot_time)
    ri0t <- psi_phi_elements$r.int.int(0,tot_time)
    Phi <- 1.0-exp(rit)/(1/f+ri0t)
    log_final_lik <- log_data_lik - log(1-Phi)
  }

  else if (cond=="crown")
  {
    rit <- psi_phi_elements$r.int(0,tot_time)
    ri0t <- psi_phi_elements$r.int.int(0,tot_time)
    Phi <- 1.0-exp(rit)/(1/f+ri0t)
    log_final_lik <- log_data_lik - log(f.lamb(tot_time)) - 2*log(1-Phi)
  }
  return(log_final_lik)
}
